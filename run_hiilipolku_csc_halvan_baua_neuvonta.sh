#!/bin/bash -l
#SBATCH --job-name=susihbaua
#SBATCH --output=/scratch/project_2002470/output/hbaua_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hbaua_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0-10# 117,217 #117,217,334 #34 #232,246,248,313,325,332,346,348 #34 #133 #34#133,134

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_baua_neuvonta.py ${SLURM_ARRAY_TASK_ID}