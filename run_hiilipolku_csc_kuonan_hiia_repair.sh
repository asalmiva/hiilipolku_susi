#!/bin/bash -l
#SBATCH --job-name=susi_khiia
#SBATCH --output=/scratch/project_2002470/output/khiia_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/khiia_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=169,269,304,369,404,469,504,569,604,669,704 #0,1,4,5,6,7,9,10,11,14,15,16,18,19,20,21,22,23,26,27,28,29,30,31,33,34,35,36,38,40,41,43,46,47,48,49,50,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,77,78,79,83,84,85,86,87,89,90,92,94,95,96,98,99

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_hiia.py ${SLURM_ARRAY_TASK_ID}