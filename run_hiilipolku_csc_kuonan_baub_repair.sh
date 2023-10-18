#!/bin/bash -l
#SBATCH --job-name=susi_kbaub
#SBATCH --output=/scratch/project_2002470/output/kbaub_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/kbaub_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=169,269,369,469,569,669 #69,304,404,504,604,704,211,311,411,511,611,711,428,528,628,728,331,431,531,631,731,635,735,337,437,537,637,737,483,583,683,391,491,591,691

#module load geoconda

module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_baub.py ${SLURM_ARRAY_TASK_ID}