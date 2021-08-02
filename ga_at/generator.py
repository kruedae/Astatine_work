

def write_sbatch(batchname):
	f = open('%s.sbatch'%batchname,'w')

	f.write("#!/bin/bash -l\n#SBATCH --nodes=1\n#SBATCH --ntasks=14\n#SBATCH --job-name=%s\n#SBATCH --partition=akanane\n#SBATCH --time=7-00:00:00\n#SBATCH --mem=140Gb\n#SBATCH --output=job.out\n#SBATCH --error=job.err\n#SBATCH --mail-user='kruedae@udel.edu'\n#SBATCH --mail-type=FAIL,END\n#SBATCH --export=NONE\nvpkg_require nwchem/6.8.1\n. /opt/shared/slurm/templates/libexec/openmpi.sh\n"%batchname
	)
	f.write('mpirun -np 18 %s.nw > %s.out'%(batchname,batchname))	
	f.close()

