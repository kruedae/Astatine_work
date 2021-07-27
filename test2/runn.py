from subprocess import Popen, PIPE
import os
import sys

def write_nw_input_file(calc_dir, temp_dir, memory, basis, pp, title, charge, mult, geometry, xcf, basename, extra_basis=None):

   ifile_name = basename + ".nw"

   f = open(ifile_name, "w")
   f.write("PERMANENT_DIR %s \n"%(calc_dir))
   f.write("SCRATCH_DIR %s \n"%(temp_dir))
   f.write("MEMORY %d mb \n"%(memory))

   f.write('\nBASIS "ao basis" spherical\n')

   if extra_basis is not None:
     f.write(extra_basis)

   for bas in basis:
      sh = bas[0]
      ex = bas[1]
      f.write(" %s   %s  \n "%("At",sh))
      f.write(" %15.7f  1.000   \n"%(ex))
   f.write("end\n")

   f.write("\n")
   f.write(pp) 

   f.write("\n")
   f.write('title  "%s" \n'%(title))
   f.write("charge %7.5f \n"%(charge))

   f.write("\n")
   f.write("geometry units angstrom\n")
   f.write("symmetry group c1\n")
   f.write(" %s \n"%(geometry))
   f.write("end\n")

   f.write("\n")
   f.write("dft \n")
   f.write("grid xfine\n")
   f.write("mult %d \n"%(mult))
   f.write("odft\n")
   f.write("iterations 100\n")
   f.write("convergence energy 1.0E-6\n")
   f.write("convergence density 5.00000E-06\n")
   f.write("convergence gradient 5.00000E-06\n")
   f.write("xc %s  \n"%(xcf))
   f.write("end \n")
   f.write("task sodft ")

   f.close()

def run_nwchem(basename):

   s = "nwchem " + basename + ".nw > " + basename + ".out"
   os.system(s)
   s = "rm " + basename + ".{movecs,db,evals}"
   os.system(s)

def get_energy(basename):

  s = "grep 'Total SO-DFT' " + basename + ".out"
  process = Popen(s, shell=True, stdout=PIPE)
  stdout, stderr = process.communicate()
  
  nwfail = False
  err = process.returncode
  if err == 0:
     line = stdout.split()
     energy = float(line[4])
  else:
     print("NWChem failed %s" % (str(err)))
     energy = -1e10
     nwfail = True

  return energy, nwfail

def run_atom(calc_dir, temp_dir, memory, basis, pp, title, charge, mult, geometry, xcf, basename):

  write_nw_input_file(calc_dir, temp_dir, memory, basis, pp, title, charge, mult, geometry, xcf, basename) 
  run_nwchem(basename) 
  
def run_mol(calc_dir, temp_dir, memory, basis, basis2, pp, title, charge, mult, geometry, xcf, basename):

  write_nw_input_file(calc_dir, temp_dir, memory, basis, pp, title, charge, mult, geometry, xcf, basename, extra_basis=basis2) 
  run_nwchem(basename) 
  
  
