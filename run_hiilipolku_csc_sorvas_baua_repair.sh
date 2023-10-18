#!/bin/bash -l
#SBATCH --job-name=susi_sbaua_test
#SBATCH --output=/scratch/project_2002470/output/sbaua_test_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/sbaua_test_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=7,8

#module load geoconda
module load python-data
# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_baua.py ${SLURM_ARRAY_TASK_ID}