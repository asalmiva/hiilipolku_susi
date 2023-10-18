#!/bin/bash -l
#SBATCH --job-name=susi_kbioa
#SBATCH --output=/scratch/project_2002470/output/kbioa_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/kbioa_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=169,194,269,294,304,369,394,404,469,494,504,569,583,594,604,669,683,684,694,704 #0,1,2,3,4,5,6,7,8,9,11,12,13,14,16,18,19,21,22,23,24,25,26,27,28,29,30,31,33,34,35,36,37,39,40,44,46,47,49,50,51,52,56,58,60,61,62,63,65,66,67,68,69,75,76,77,79,80,81,83,84,88,90,91,94,95,99

#module load geoconda

module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_bioa.py ${SLURM_ARRAY_TASK_ID}