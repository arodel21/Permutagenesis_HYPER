#!/bin/bash
#SBATCH --partition=<node_name>
#SBATCH --gpus=<n_gpus>
#SBATCH --job-name=<job_name>
#SBATCH --ntasks=<n_threads>
#SBATCH --mem=50G
#SBATCH --time=00:10:00
#SBATCH --mail-user=<user_email>
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Define output directory based on SLURM variables.
OUT_DIR="Sbatches/${SLURM_JOB_NAME}"
mkdir -p "$OUT_DIR"

# Redirect stdout and stderr to files in the new folder.
exec > "$OUT_DIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
exec 2> "$OUT_DIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"

module load python-cbrg
module load bedtools

#JSON=07_JSONs/Human_WG/chr1\:1-196609.json
JSON=$1

python 00_Resources/pipeline.py --json_file $JSON
