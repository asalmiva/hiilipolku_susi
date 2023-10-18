#!/bin/bash -l
#SBATCH --job-name=susikbaua
#SBATCH --output=/scratch/project_2002470/output/k_baua_job_out_%A_%a.txt
#SBATCH --error=/scratch/project_2002470/error/k_baua_job_err_%A_%a.txt
#SBATCH --account=project_2002470
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000
#SBATCH --array=20,24,44,48,33,84,8,74,30,17,9,40,87,41,99,92,70,90,75,151,150,153,122,146,182,100,135,178,191,196,113,163,159,107,104,106,227,221,223,243,283,298,289,257,276,285,256,258,261,124,225,336,347,187,352,328,332,381,377,329,386,302,175,190,244,426,251,117,233,427,493,466,200,412,320,436,263,405,422,321,396,443,485,516,515,562,447,642,639,597,501,584,413,474,409,428,499,224,378,695,671,550,648,557,300,586,541,421,556,692,676,630,420,640,506,544,633,490,651,575,621,500,733,70,191,178,135,104,283,381,378,681,691,570,469,691,681,570,469,469,570,691,469 #0-99

#module load geoconda
module load python-data
# pip install --user openpyxl # tarvitaanko tämä?

python /scratch/project_2002470/SUSI_HIILIPOLKU/run_hiilipolku_csc_kuonan_baua.py ${SLURM_ARRAY_TASK_ID}