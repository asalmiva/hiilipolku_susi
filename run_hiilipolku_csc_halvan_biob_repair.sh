#!/bin/bash -l
#SBATCH --job-name=susi_baub
#SBATCH --output=/scratch/project_2002470/output/baub_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/baub_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0,1,4,8,11,13,21,22,23,24,25,27,28,29,30,32,35,36,39,41,43,46,48,51,52,53,54,55,57,58,59,61,62,63,64,66,67,68,77,78,85,86,95,96

module load geoconda

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob.py ${SLURM_ARRAY_TASK_ID}