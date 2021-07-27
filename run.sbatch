#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --job-name=GA_DZ_PP
#SBATCH --partition=akanane
#SBATCH --time=7-00:00:00
#SBATCH --mem=140Gb
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --mail-user='akanane@udel.edu'
#SBATCH --mail-type=FAIL,END
#SBATCH --export=NONE

vpkg_require nwchem/6.8.1
. /opt/shared/slurm/templates/libexec/openmpi.sh

python3 main.py 

