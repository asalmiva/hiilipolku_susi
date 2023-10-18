#!/bin/bash -l
#SBATCH --job-name=susi_baub
#SBATCH --output=/scratch/project_2002470/output/baub_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/baub_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0,1,9,11,23,24,25,26,27,28,30,32,33,36,39,44,46,47,48,49,52,53,54,55,57,58,59,62,63,66,70,71,78,79,82,83,86,92,95

module load geoconda

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_bioa.py ${SLURM_ARRAY_TASK_ID}