#!/bin/bash
#$ -l h_vmem=40G
#$ -l mem_free=40G
#$ -t 1-137
#$ -l h_rt=48:00:00
#$ -pe smp 1
#$ -R yes
#$ -V

#________________________________________SOURCE PHENIX/QFIT_____________________________________________#
source /wynton/home/fraserlab/jessicaflowers/phenix/phenix-1.20.1-4487/phenix_env.sh
export PATH="/wynton/home/fraserlab/jessicaflowers/miniconda3/bin:$PATH"

source activate qfit_ligand
which python
export PHENIX_OVERWRITE_ALL=true

#________________________________________________PDB INFO________________________________________________#

PDB_file=/wynton/home/fraserlab/jessicaflowers/tp_qfit_ligands/tp_ligands.txt
PDB_dir='/wynton/home/fraserlab/jessicaflowers/tp_qfit_ligands'
output_dir='/wynton/home/fraserlab/jessicaflowers/tp_qfit_ligands/output_data/results_from_qfit_run_1'

PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)
category='qfit_run_1'

# Extract and format the required information
pdb_id=$(echo "$PDB" | awk -F, '{print $1}')
path=${PDB_dir}/${pdb_id}/${category}
chain_res=$(echo "$PDB" | awk -F, '{print $(NF-3) "," $(NF-2)}')

cd ${path}


#________________________________________________SPLIT MULTICONF________________________________________________
echo ${pdb_id} 
echo ${path}

# Split the refined qFit-Ligand output multiconformer model into seperate PDBs

# this should save multi_ligand_A.pdb, multi_ligand_B.pdb, ...
qfit_split_conf=$(split_multiconformer_ligand.py ${pdb_id}_qFit_ligand.pdb  --residue=${chain_res} --directory=${path} --output_name qfit)

# this should save depo_ligand_A.pdb and depo_ligand_B.pdb
depo_split_conf=$(split_multiconformer_ligand.py ${pdb_id}.pdb  --residue=${chain_res} --directory=${path} --output_name depo)


#________________________________________________METRICS________________________________________________
cd ${output_dir}

# RSCC

if [[ -f "${path}/${pdb_id}_qFit_ligand.pdb" ]]; then
    rscc=$(compare_rscc_voxel.py ${path}/${pdb_id}_001.pdb ${path}/composite_omit_map.mtz --gen_pdb ${path}/${pdb_id}_qFit_ligand.pdb --gen_map ${path}/${pdb_id}_qFit_ligand.mtz --residue ${chain_res} --pdb ${pdb_id} --directory ${output_dir})
else
    rscc=$(compare_rscc_voxel.py ${path}/${pdb_id}_001.pdb ${path}/composite_omit_map.mtz --gen_pdb ${path}/multiconformer_ligand_bound_with_protein.pdb --gen_map ${path}/composite_omit_map.mtz --residue ${chain_res} --pdb ${pdb_id} --directory ${output_dir})
fi

# Parse the log file for total run time, and number of output conformers
python ${PDB_dir}/parse_log_file.py --pdb ${pdb_id} --path ${path} --output_dir ${output_dir}

# RMSD between qfit output conformers 
python ${PDB_dir}/calc_rmsd.py --pdb ${pdb_id} --path ${path} --output_dir ${output_dir} --conf_type qfit_ligand

# Torsion strain of qfit conformers and deposited conformers


python ${PDB_dir}/calc_torsion_strain.py --pdb ${pdb_id}_qfit --path ${path} --output_dir ${output_dir} --conf_type qfit_ligand
python ${PDB_dir}/calc_torsion_strain.py --pdb ${pdb_id}_depo --path ${path} --output_dir ${output_dir} --conf_type depo_ligand

