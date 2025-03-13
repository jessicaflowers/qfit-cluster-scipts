#!/bin/bash
#$ -l h_vmem=4G
#$ -l mem_free=4G
#$ -t 1-100 # change based on how many strucutres in dataset 
#$ -l h_rt=28:00:00
#$ -pe smp 8

# Pre-QFit Refinement Script

#________________________________________________INPUTS________________________________________________#
PDB_file="/wynton/home/fraserlab/jessicaflowers/survey_true_pos/tp_inputs.txt"
base_dir="/wynton/home/fraserlab/jessicaflowers/survey_true_pos"

#________________________________________________SET PATHS________________________________________________#
source /wynton/home/fraserlab/jessicaflowers/phenix/phenix-1.20.1-4487/phenix_env.sh
export PATH="/wynton/home/fraserlab/jessicaflowers/miniconda3/bin:$PATH"
source activate lig_dev
export PHENIX_OVERWRITE_ALL=true

#________________________________________________UTILITY FUNCTIONS________________________________________________#
log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $*"
}

run_command() {
    log "Running: $*"
    eval "$*"
    local status=$?
    if [ $status -ne 0 ]; then
        log "Error: Command failed - $*"
        exit $status
    fi
}

# Get xray_data_labels by checking available field names in the MTZ file.
get_xray_data_labels() {
    local mtz_file="$1"
    local mtzmetadata
    mtzmetadata=$(phenix.mtz.dump "${mtz_file}")
    local obstypes=("FP" "FOBS" "F-obs" "I" "IOBS" "I-obs" "F(+)" "I(+)" "FSIM")
    local label=""
    for field in "${obstypes[@]}"; do
        if grep -F -q -w "SIG${field}" <<< "${mtzmetadata}"; then
            label="${field},SIG${field}"
            break
        fi
    done
    if [ -z "$label" ]; then
        echo >&2 "Could not determine Fo field name with corresponding SIGFo in ${mtz_file}."
        echo >&2 "Checked fields: ${obstypes[*]}"
        exit 1
    fi
    echo "$label"
}

# Run phenix.refine using the provided parameters.
run_refinement() {
    local pdb_id="$1"
    local pdb_input="$2"
    local cif_input="$3"
    local nqh_flip="$4"
    local labels="$5"
    phenix.refine "${pdb_id}.mtz" "${pdb_input}" \
      refinement.input.monomers.file_name="${cif_input}" \
      refinement.refine.strategy=individual_sites+individual_adp+occupancies \
      refinement.output.prefix="${pdb_id}" \
      refinement.main.number_of_macro_cycles=5 \
      refinement.main.nqh_flips="${nqh_flip}" \
      refinement.output.write_maps=False \
      refinement.hydrogens.refine=riding \
      refinement.main.ordered_solvent=True \
      refinement.target_weights.optimize_xyz_weight=true \
      refinement.target_weights.optimize_adp_weight=true \
      refinement.input.xray_data.r_free_flags.generate=True \
      refinement.input.xray_data.labels="${labels}"
}

#________________________________________________RUN PRE-QFIT REFINEMENT________________________________________________#
# Get the current PDB entry using sed (more concise than cat/head/tail)
PDB=$(sed -n "${SGE_TASK_ID}p" "${PDB_file}")
pdb_id=$(echo "$PDB" | awk -F, '{print $1}')
lig_name=$(echo "$PDB" | awk -F, '{print $NF}')

log "Ligand name: ${lig_name}"

cd "${base_dir}/${pdb_id}" || { log "Directory not found: ${base_dir}/${pdb_id}"; exit 1; }
log "Working in directory ${base_dir}/${pdb_id}"

# Remove alternative conformations and duplicates
run_command "remove_altconfs ${pdb_id}.pdb"
run_command "remove_duplicates ${pdb_id}.single.pdb"

# Convert CIF to MTZ and move the file
run_command "phenix.cif_as_mtz ${pdb_id}-sf.cif --merge"
if [ ! -f "${pdb_id}-sf.mtz" ]; then
    log "Merge failed for ${pdb_id}-sf.cif"
fi
log "Moving ${pdb_id}-sf.mtz to ${pdb_id}.mtz"
mv "${pdb_id}-sf.mtz" "${pdb_id}.mtz"

# Determine the X-ray data labels from the MTZ file
xray_data_labels=$(get_xray_data_labels "${pdb_id}.mtz")
log "Data labels: ${xray_data_labels}"

# Run ready_set to generate a fixed PDB file
log "Running ready_set"
run_command "phenix.ready_set ${pdb_id}.single.pdb.fixed"

# Run refinement based on the existence of a .ligands.cif file, else try elbow
if [ -f "${pdb_id}.single.pdb.ligands.cif" ]; then
    log "Using ligands.cif for refinement"
    run_refinement "${pdb_id}" "${pdb_id}.single.pdb.updated.pdb" "${pdb_id}.single.pdb.ligands.cif" "True" "${xray_data_labels}"
else
    log "ligands.cif not found; trying phenix.elbow"
    run_command "phenix.elbow ${pdb_id}.single.pdb.fixed --residue ${lig_name}"
    elbow_cif="elbow.${lig_name}.${pdb_id}_single_pdb_fixed.cif"
    if [ -f "${elbow_cif}" ]; then
        run_refinement "${pdb_id}" "${pdb_id}.single.pdb.fixed" "${elbow_cif}" "True" "${xray_data_labels}"
    else
        log "Error: phenix.elbow failed to generate ${elbow_cif}"
    fi
fi

# If a reduce error is detected, re-run refinement with adjusted parameters
if [ -f "reduce_failure.pdb" ]; then
    if [ -f "${pdb_id}.single.pdb.ligands.cif" ]; then
        log "Reduce failure detected; re-running refinement with ligands.cif (nqh_flips=False)"
        run_refinement "${pdb_id}" "${pdb_id}.single.pdb.updated.pdb" "${pdb_id}.single.pdb.ligands.cif" "False" "${xray_data_labels}"
    else
        log "Reduce failure detected; trying elbow with nqh_flips=False"
        run_command "phenix.elbow ${pdb_id}.single.pdb.fixed --residue ${lig_name}"
        elbow_cif="elbow.${lig_name}.${pdb_id}_single_pdb_fixed.cif"
        if [ -f "${elbow_cif}" ]; then
            run_refinement "${pdb_id}" "${pdb_id}.single.pdb.fixed" "${elbow_cif}" "False" "${xray_data_labels}"
        else
            log "Error: phenix.elbow failed in reduce_failure block"
        fi
    fi
fi

# Run composite omit map
log "Running composite omit map"
phenix.composite_omit_map "${pdb_id}.mtz" "${pdb_id}_001.pdb" omit-type=refine nproc="${NSLOTS:-1}" r_free_flags.generate=True exclude_bulk_solvent=True
