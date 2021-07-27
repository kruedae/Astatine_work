import os
import sys
import multiprocessing
from multiprocessing import Process
import numpy as np
import random

import basis
import runn
import algo

Ha_to_eV = 27.2114

n_gen = 5000
n_ind = 100

save_n_best = 10.0  # % of the entire generation will be saved from the tournament

fittest_win = 0.8   # is the probability that better performing basis wins the face-to-face comp.

job_type = 0
save_file = "job_save.txt"

#------------ Reference data ---------------
E_At  = -261.958026856090    # SO HF/uncontracted-cc-pV5Z-PP atomic calc. starting from anion orbitals
E_Atp = -261.653246827283    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
E_Atn = -262.007008940345    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
E_At2 = -526.187655911062 
E_HAt = -262.502879716577
E_AtI = -559.678938556197

#At_EA = 2.412
#At_IP = 9.31751

#-------------------------------------------

input_basis_file = "/work/akanane/users/akanane/Relativistic/GA/cc-pVDZ-F-pp_uncontracted.bas"
H_basis_file     = "/work/akanane/users/akanane/Relativistic/GA/H_uccpv5z.bas"
I_basis_file     = "/work/akanane/users/akanane/Relativistic/GA/I_uccpv5z.bas"
AtI_pp_file      = "/work/akanane/users/akanane/Relativistic/GA/AtI_gatchina.pp"
At_pp_file       = "/work/akanane/users/akanane/Relativistic/GA/At_gatchina.pp"

cal_dir          = os.getcwd()
scratch_dir1      = "/lustre/scratch/akanane/nwchem/1"
scratch_dir2      = "/lustre/scratch/akanane/nwchem/2"
scratch_dir3      = "/lustre/scratch/akanane/nwchem/3"
scratch_dir4      = "/lustre/scratch/akanane/nwchem/4"
scratch_dir5      = "/lustre/scratch/akanane/nwchem/5"
scratch_dir6      = "/lustre/scratch/akanane/nwchem/6"

memory     = 32000   # in mb
xcf        = "pbe0"
xhf        = "HFexch"

At_title   = "At"
Atn_title  = "At-"
Atp_title  = "At+"
At2_title  = "At2"
HAt_title  = "HAt"
AtI_title  = "AtI"

At_charge  = 0
Atn_charge = -1
Atp_charge = 1
At2_charge = 0
HAt_charge = 0
AtI_charge = 0

At_mult    = 2
Atn_mult   = 1
Atp_mult   = 3
At2_mult   = 1
HAt_mult   = 1
AtI_mult   = 1

At_geometry = "At   0.0   0.0   0.0"
At2_geometry= "At  -1.48753783     0.00000000     0.00000000\n At    1.48753783     0.00000000     0.00000000"
HAt_geometry= "H   -1.69812258     0.00000000     0.00000000\n At    0.01145529     0.00000000     0.00000000"
AtI_geometry= "At  -1.12002211     0.00000000     0.00000000\n I     1.69973226     0.00000000     0.00000000" 

basenameAt  = "nwc"
basenameAtp = "nwcp"
basenameAtn = "nwcn"
basenameAt2 = "nwcd"
basenameHAt = "nwchat"
basenameAtI = "nwcati"

emin = 1e20
gens = -1
inds = -1
eoriginal = 0.0
#
# read input basis to be modified and print it out
#

if __name__=="__main__":

   if not os.path.exists(scratch_dir1):
      os.mkdir(scratch_dir1)
   if not os.path.exists(scratch_dir2):
      os.mkdir(scratch_dir2)
   if not os.path.exists(scratch_dir3):
      os.mkdir(scratch_dir3)
   if not os.path.exists(scratch_dir4):
      os.mkdir(scratch_dir4)
   if not os.path.exists(scratch_dir5):
      os.mkdir(scratch_dir5)
   if not os.path.exists(scratch_dir6):
      os.mkdir(scratch_dir6)

   out_file_name = "ga.out"
   fo = open(out_file_name, "w")

   if job_type == 0:
      fo.write(" Starting new calculation ")
      bas_start = basis.read_start_basis(input_basis_file)

      fo.write(" -------- Input basis -------- \n")
      for p in range(len(bas_start)):
         v = bas_start[p]
         fo.write("   %s     %18.9f \n"%(v[0], v[1]))
      fo.write(" ------ Input basis end ------ \n\n")

   #
   # read pseudopotential
   #
   fo.write(" Reading pseudopotential from %s "%(At_pp_file))
   At_pp = basis.read_pseudopotential(At_pp_file)

   fo.write(" -------- Pseudopotential ----------\n")
   fo.write(At_pp)
   fo.write(" -------- Pseudopotential end ----------\n\n")

   #
   # read H and I basis sets and pseudopotential
   # 
   H_basis = basis.read_atom_basis(H_basis_file)
   I_basis = basis.read_atom_basis(I_basis_file)
   AtI_pp  = basis.read_pseudopotential(AtI_pp_file)  

   #
   # clone input basis into n_ind 
   #
   original_basis_mut = True
   basis_list  = algo.clone(bas_start, n_ind)
   bas_mutated, idxn = algo.mutate(basis_list, original_basis_mut)

   for gen in range(n_gen):

      if gen == 0 and job_type == 1:
         fo.write(" Restarting previous job. Reading basis sets from \n"%(save_file))
         bas_mutated = basis.read('job_saved.txt')


      fo.write("\n\n <<<<<<<<<<<<<<<< Beginning of a new generation %d out of %d >>>>>>>>>>>>>>>>>>>\n\n"%((gen+1),n_gen))
      fo.flush()

      n_failed = 0
      n_cur_gen = len(bas_mutated)
      fo.write(" Current generation size %d \n"%(n_cur_gen))

      errors = np.zeros((n_cur_gen))
      ips = np.zeros((n_cur_gen))
      eas = np.zeros((n_cur_gen))

      for n, bas in enumerate(bas_mutated):

         # 
         # run NWChem calculations for all individuals in this generation
         # 
         # At atom
         p1 = Process(target=runn.run_atom, args=(cal_dir, scratch_dir1, memory, bas, 
                                                  At_pp, At_title, At_charge, At_mult, 
                                                  At_geometry, xhf, basenameAt,))
         # At+
         p2 = Process(target=runn.run_atom, args=(cal_dir, scratch_dir2, memory, bas, 
                                                  At_pp, Atp_title, Atp_charge, Atp_mult, 
                                                  At_geometry, xhf, basenameAtp,))
         # At-
         p3 = Process(target=runn.run_atom, args=(cal_dir, scratch_dir3, memory, bas, 
                                                  At_pp, Atn_title, Atn_charge, Atn_mult, 
                                                  At_geometry, xhf, basenameAtn,))

         # At2
         #p4 = Process(target=runn.run_atom, args=(cal_dir, scratch_dir4, memory, bas,
         #                                         At_pp, At2_title, At2_charge, At2_mult,
         #                                         At2_geometry, xcf, basenameAt2,))
         # HAt
         #p5 = Process(target=runn.run_mol, args=(cal_dir, scratch_dir5, memory, bas,
         #                                        H_basis, At_pp, HAt_title, HAt_charge, 
         #                                        HAt_mult, HAt_geometry, xcf, 
         #                                        basenameHAt,))

         # IAt
         #p6 = Process(target=runn.run_mol, args=(cal_dir, scratch_dir6, memory, bas,
         #                                        I_basis, AtI_pp, AtI_title, AtI_charge, 
         #                                        AtI_mult, AtI_geometry, xcf, 
         #                                        basenameAtI,))

         p1.start()
         p2.start()
         p3.start()
         #p4.start()
         #p5.start()
         #p6.start()
         
         p1.join()
         p2.join()
         p3.join()
         #p4.join()
         #p5.join()
         #p6.join()

         E_out_at,  nwfail1 = runn.get_energy(basenameAt)
         E_out_atp, nwfail2 = runn.get_energy(basenameAtp)
         E_out_atn, nwfail3 = runn.get_energy(basenameAtn)
         #E_out_at2, nwfail4 = runn.get_energy(basenameAt2)
         #E_out_hat, nwfail5 = runn.get_energy(basenameHAt)
         #E_out_ati, nwfail6 = runn.get_energy(basenameAtI)

         #
         #--------------------------------------------------------------------------- 

         if nwfail1 or nwfail2 or nwfail3: #or nwfail4 or nwfail5 or nwfail6:
            n_failed += 1 

         errors[n]  = np.abs(E_At  - E_out_at) 
         errors[n] += np.abs(E_Atp - E_out_atp) 
         errors[n] += np.abs(E_Atn - E_out_atn)
         #errors[n] += np.abs(E_At2 - E_out_at2)
         #errors[n] += np.abs(E_HAt - E_out_hat)
         #errors[n] += np.abs(E_AtI - E_out_ati)    

         errors[n] *= Ha_to_eV/3.0 

         #ips[n] = Ha_to_eV*(E_out_atp - E_out_at)
         #eas[n] = Ha_to_eV*(E_out_at - E_out_atn)

 
         err_diff = np.max(errors) - np.min(errors)

         #
         # first time going through this let's highlight the original basis 
         #
         if n==idxn and original_basis_mut:
            fo.write(" ** idx = %d error = %10.7f [eV] *** this is the original basis. \n"%(n,errors[n]))
            original_basis_mut = False
            eoriginal = errors[n]
         else:
            fo.write(" ** idx = %d error = %10.7f [eV] \n"%(n,errors[n]))
         fo.flush()
   
      fo.write("\n\n Done running NWChem calculations. \n\n")
      fo.flush()
   
      #
      # sort all results by everage error
      #
      ind = np.argsort(errors)

      #
      # let's deal with failed calculations first
      #  
      if n_failed > 0:
         fo.write(" NWChem failed for the following basis sets (possibly due to linear dependences): \n")
         for nx in range(n_failed):
            fo.write(" %d "%(ind[n_cur_gen-1-nx]))
         fo.write("\n")
      fo.write(" Total calculations failed : %d \n"%(n_failed))
 
      #
      # do some accounting here...
      # 
      n_rem     = n_cur_gen - n_failed
      n_promote = int(round(n_rem*1.0/save_n_best))
      n_tourn   = n_rem - n_promote

      if n_tourn % 2 == 1:
         n_promote += 1
         n_tourn   -= 1     # will make it even

      n_cross = n_ind - n_promote - n_tourn // 2
   
      fo.write(" Directly promoted to the next generation : %d \n"%(n_promote))
      fo.write(" Going into the tournament stage : %d \n"%(n_tourn))

      fo.write(" Best performing basis in this generation gives the error of %10.7f [eV] \n"%(errors[ind[0]]))
      #fo.write(" IP = %7.5f [eV], EA = %7.5f [eV] \n"%(ips[ind[0]], eas[ind[0]]))
      #fo.write(" Error in IP w.r.t. experiment = %7.5f [eV] \n"%(ips[ind[0]] - At_IP))
      #fo.write(" Error in EA w.r.t. PRA 91, 020501 (2015) = %7.5f [eV] \n"%(eas[ind[0]] - At_EA))
      basis.write_basis(bas_mutated[ind[0]], gen, ind[0])
      fo.write(" See file bas_gen%d_id%d.dat \n"%(gen,ind[0]))

      fo.write(" Difference between best and worst performing bases at %d generation is: %10.7f \n"%(gen,err_diff))

      if errors[ind[0]] < emin:
         de = errors[ind[0]] - emin
         emin = errors[ind[0]]
         gens = gen
         inds = ind[0]
         if gen > 0:
           fo.write(" Smallest MAE updated: gen = %d value = %7.5f [eV], diff = %7.5f [eV], file = bas_gen%d_id%d.dat \n"%(gen,emin,de,gens,inds))
           fo.write(" Difference w.r.t. original basis = %10.7f [eV] \n"%(eoriginal-emin))
      else:
         fo.write(" Smallest error has not improved since previous step: %10.7f [eV] \n"%(emin))

      #
      #  print a few basis sets giving the lowest errors in the current generation
      #
      fo.write(" %d lowest error indices are: \n"%(n_promote))
      fo.flush()
      for nx in range(n_promote):
        fo.write(" %d "%(ind[nx]))
      fo.write("\n")

      #
      # -- Tournament
      #
      fo.write("\n\n %d basis sets are going into the tournament: "%(n_tourn))

      ind_tourn = ind[n_promote:n_cur_gen-n_failed]
      for nx in range(n_tourn):
         fo.write(" %d "%(ind_tourn[nx]))
      fo.write("\n\n")
      fo.flush()

      #
      # save n best basis sets
      #
      next_gen = []
      for nx in range(n_promote):
         idd = ind[nx] 
         next_gen.append(bas_mutated[idd])

      # tournament stage:
      np.random.shuffle(ind_tourn)
      n_pairs = n_tourn // 2
      fo.write(" Tournament begins...there will be %d pairs \n"%(n_pairs))
      fo.flush()
      for rnd in range(n_pairs):
         id1 = ind_tourn[2*rnd]
         id2 = ind_tourn[2*rnd+1]
         er1 = errors[id1]
         er2 = errors[id2]
         rg = np.random.random()
         fo.write(" Tournament %d between %d and %d with errors %7.5f and %7.5f, random = %7.5f ..."%(rnd,id1,id2,er1,er2,rg)) 
         if er1 < er2:
            if rg < fittest_win:
               winner = id1
            else:
               winner = id2
         else:
            if rg < fittest_win:
               winner = id2
            else:
               winner = id1
         fo.write(" winner = %d \n"%(winner)) 
         next_gen.append(bas_mutated[winner])
      fo.write(" Tournament ends.\n\n")
      fo.flush()

      # now we will do the mating to restore the original size of the population
      fo.write(" we will now start crossover stage to add %d more individuals to current generation.\n"%(n_cross))
      fo.flush()

      ind_here = list(range(n_cross))
      random.shuffle(ind_here)
      n_mat = n_cross // 2
      for nx in range(n_mat):
        ind1 = ind_here[2*nx]
        ind2 = ind_here[2*nx+1]
        mat_out1, mat_out2 = algo.crossover(next_gen[ind1], next_gen[ind2])
        next_gen.append(mat_out1)
        next_gen.append(mat_out2)
        fo.write(" crossover between %d and %d \n"%(ind1,ind2))
   
      # mutate and go back to the top of the loop
      fo.write(" mutating ....\n")
      bas_mutated, idxn = algo.mutate(next_gen, False)   
   
      # save checkpoint
      fo.write(" Saving active basis sets to %s \n"%(save_file))
      basis.save(bas_mutated, save_file)

   fo.close()
