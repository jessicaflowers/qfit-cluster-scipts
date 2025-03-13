#!/bin/bash
#$ -l h_vmem=40G
#$ -l mem_free=40G
#$ -t 1-133
#$ -l h_rt=48:00:00
#$ -pe smp 1
#$ -R yes
#$ -V

#________________________________________SOURCE PHENIX/QFIT_____________________________________________#
source /wynton/home/fraserlab/jessicaflowers/phenix/phenix-1.20.1-4487/phenix_env.sh
export PATH="/wynton/home/fraserlab/jessicaflowers/miniconda3/bin:$PATH"

source activate qfit_refine
which python
export PHENIX_OVERWRITE_ALL=true

#________________________________________FILE AND DIRECTORY PATHS_____________________________________________#
PDB_file="/wynton/home/fraserlab/jessicaflowers/high_strain_evaluation/strain_evaluation_with_ligs_clean.txt"
PDB_dir="/wynton/home/fraserlab/jessicaflowers/high_strain_evaluation"
OUTPUT_DIR="/wynton/home/fraserlab/jessicaflowers/high_strain_evaluation/output_data/unref_kcal_mol_10.txt"
CALC_SCRIPT="/wynton/home/fraserlab/jessicaflowers/high_strain_evaluation/calc_torsion_strain.py"

#________________________________________LOOP OVER EACH LINE_____________________________________________#
while IFS=',' read -r pdb chain resnum lig; do
  # Remove any extra whitespace
  pdb=$(echo "$pdb" | xargs)
  chain=$(echo "$chain" | xargs)
  resnum=$(echo "$resnum" | xargs)
  lig=$(echo "$lig" | xargs)

  # Construct the ligand file identifier 
  conf_type="${pdb}_${lig}_${chain}_${resnum}"

  # Determine the directory where the PDB file is located.
  path="${PDB_dir}/${pdb}"
  echo "Processing ligand ${conf_type} in directory ${path}"

  # Change to the pdb folder (which contains the ligand pdb file)
  cd "$path" || { echo "Directory ${path} not found. Skipping ${conf_type}."; continue; }

  # Run the strain calculation.
  # This command expects that there is a ligand PDB file named "${conf_type}.pdb" in the current directory.
  python "$CALC_SCRIPT" --pdb "$pdb" --path "$path" --output_dir "$OUTPUT_DIR" --conf_type "$conf_type"

done < "$PDB_file"

echo "All strain calculations complete."
