#!/bin/bash -l
#SBATCH --job-name=susishiia
#SBATCH --output=/scratch/project_2002470/output/shiia_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/shiia_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0-99

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?pip install --user XlsxWriter

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_sorvas_hiia.py ${SLURM_ARRAY_TASK_ID}