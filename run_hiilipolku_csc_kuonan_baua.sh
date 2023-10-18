#!/bin/bash -l
#SBATCH --job-name=susikbaua
#SBATCH --output=/scratch/project_2002470/output/k_baua_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/k_baua_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0-99

#module load geoconda
module load python-data
# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_baua.py ${SLURM_ARRAY_TASK_ID}