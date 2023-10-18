#!/bin/bash -l
#SBATCH --job-name=susi_hiib
#SBATCH --output=/scratch/project_2002470/output/hiib_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hiibb_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0,1,4,5,9,10,13,15,21,23,25,27,28,30,32,34,41,43,44,45,46,47,48,50,55,56,57,58,59,60,61,62,63,64,65,66,68,69,71,76,77,82,83,85,92,95,96,98

module load geoconda

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_hiib.py ${SLURM_ARRAY_TASK_ID}