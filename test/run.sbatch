#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=GA_TZ_PP
#SBATCH --partition=akanane
#SBATCH --time=7-00:00:00
#SBATCH --mem=140Gb
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --mail-user='akanane@udel.edu'
#SBATCH --mail-type=FAIL,END
#SBATCH --export=NONE

vpkg_require nwchem
python main.py 

