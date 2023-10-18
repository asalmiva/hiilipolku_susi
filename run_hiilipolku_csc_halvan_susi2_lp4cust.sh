#!/bin/bash -l
#SBATCH --job-name=susihbaua
#SBATCH --output=/scratch/project_2002470/output/hbaua_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hbaua_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=5 # 0-10# 117,217 #117,217,334 #34 #232,246,248,313,325,332,346,348 #34 #133 #34#133,134

module load python-data

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_lp4_cust.py ${SLURM_ARRAY_TASK_ID}
#python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_neuvonta.py ${SLURM_ARRAY_TASK_ID}
#python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_lp4.py ${SLURM_ARRAY_TASK_ID}
#python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_vol01.py ${SLURM_ARRAY_TASK_ID}
#python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_vol001.py ${SLURM_ARRAY_TASK_ID}
#python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_biob_volnage.py ${SLURM_ARRAY_TASK_ID}
