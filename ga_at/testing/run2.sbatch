#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks=60
#SBATCH --job-name=nodeprobe
#SBATCH --partition=akanane
#SBATCH --time=7-00:00:00
#SBATCH --mem=140Gb
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --mail-user='kruedae@udel.edu'
#SBATCH --mail-type=FAIL,END
#SBATCH --export=NONE

vpkg_require nwchem/6.8.1
. /opt/shared/slurm/templates/libexec/openmpi.sh

mpirun -np 60 nwchem nwcat2.nw > testing60cores.out


