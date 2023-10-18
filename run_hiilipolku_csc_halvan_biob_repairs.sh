#!/bin/bash -l
#SBATCH --job-name=susihbiob
#SBATCH --output=/scratch/project_2002470/output/hbiob_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hbiob_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=35 #130 #,217 #35 #153,232,246,248,253,325,332,346,348,353 #90 #317 #180,217,315,316 #217 #91 #6,20,22,40,60,65,73,74,94,110,114,117,142,151,164,169,176,180,190,196,197,198,207,213,217,218,222,267,275,280,293,305,306,307,315,316,342,351,356,365,367,369 #40 #133,134

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob.py ${SLURM_ARRAY_TASK_ID}