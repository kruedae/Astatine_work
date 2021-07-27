#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --job-name=GA_TZ_n
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

python3 /work/akanane/users/akanane/Relativistic/GA/parallel/main.py input.json

