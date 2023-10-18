#!/bin/bash -l
#SBATCH --job-name=susihhiib
#SBATCH --output=/scratch/project_2002470/output/hhiib_array_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/hhiib_array_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=100,101,113,121,125,132,156,157,159,162,163,164,168,169,182,185,192,196,200,201,204,205,213,215,221,223,225,228,232,244,246,248,256,257,259,260,261,262,263,264,265,268,269,271,276,277,282,283,285,292,295,296,298,300,301,304,305,309,310,313,315,321,323,325,327,328,330,332,341,343,344,345,346,347,348,350,355,356,357,358,359,360,361,362,363,364,365,366,368,369 #34 #134

#module load geoconda
module load python-data

# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_halvan_hiib.py ${SLURM_ARRAY_TASK_ID}