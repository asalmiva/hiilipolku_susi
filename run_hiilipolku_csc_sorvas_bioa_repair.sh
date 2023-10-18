#!/bin/bash -l
#SBATCH --job-name=susi_baub
#SBATCH --output=/scratch/project_2002470/output/baub_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/baub_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0-99

module load geoconda

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_bioa.py ${SLURM_ARRAY_TASK_ID}