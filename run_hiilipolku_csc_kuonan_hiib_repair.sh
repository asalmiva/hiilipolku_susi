#!/bin/bash -l
#SBATCH --job-name=susi_baub
#SBATCH --output=/scratch/project_2002470/output/baub_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/baub_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=169,194,269,294,304,369,394,399,404,469,494,499,504,569,583,594,599,604,669,683,684,694,699,704 #0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,33,34,35,36,37,38,39,40,41,43,44,45,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,77,78,79,81,83,84,85,87,88,91,92,94,96,98,99

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_hiib.py ${SLURM_ARRAY_TASK_ID}