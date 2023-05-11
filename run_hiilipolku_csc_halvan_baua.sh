#!/bin/bash -l
#SBATCH --job-name=susi
#SBATCH --output=/scratch/project_2002470/output/array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-99

module load geoconda

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_baua.py ${SLURM_ARRAY_TASK_ID}