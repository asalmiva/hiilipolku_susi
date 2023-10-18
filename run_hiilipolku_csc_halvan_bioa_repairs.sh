#!/bin/bash -l
#SBATCH --job-name=susihbioa
#SBATCH --output=/scratch/project_2002470/output/hbioa_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hbioa_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=217 #232,246,248,325,332,346,348 #90 #317 #180,217,315,316 #217 #91 #6,20,22,40,60,65,73,74,94,110,114,117,142,151,164,169,176,180,190,196,197,198,207,213,217,218,222,267,275,280,293,305,306,307,315,316,342,351,356,365,367,369 #40 #133,134

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_bioa.py ${SLURM_ARRAY_TASK_ID}