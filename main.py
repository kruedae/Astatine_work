import os
import sys
import multiprocessing
from multiprocessing import Process
import numpy as np
import random

import basis
import runn
import algo

import argparse
import json

def main():

   #------------ Reference data ----------
   Ha_to_eV = 27.211
   E_At  = -263.0812543    # SO HF/uncontracted-cc-pV5Z-PP atomic calc. starting from anion orbitals
   E_Atp = -262.7419035    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
   E_At2p = -262.0950157    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
   E_Atn = -263.1663307    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
   E_At_m_calc = np.array([E_At2p, E_Atp, E_At, E_Atn])
   # Put the values of distance and energy
   R_At2 = [2.97842834, 2.77842834, 2.87842834, 3.07842834, 3.17842834 ]
   E_At2 = np.array([-526.1925101, -526.1867749, -526.191208, -526.1916841, -526.1894469])
   R_HAt = [1.5237857, 1.6237857, 1.7237857, 1.8237857, 1.9237857]
   E_HAt = np.array([-263.6547155, -263.6666387, -263.670009, -263.6677744, -263.6619657])
   N_At = len(E_At_m_calc)
   N_At2 = len(E_At2)
   N_HAt = len(E_HAt)
   #At_EA = 2.412
   #At_IP = 9.31751
   # Use parser function to write in console (review)
   parser = argparse.ArgumentParser()
   parser.add_argument("input_file",   help="input json file")
   args = parser.parse_args()

   # read json file
   with open(args.input_file) as json_file:
      j = json.load(json_file)

   # read parameters from json file
   n_gen            = j['generations']
   n_ind            = j['population_size']
   save_n_best      = j['save_n_best']
   fittest_win      = j['fittest_win_probability']/100.0
   job_type         = j['job_type']
   save_file        = j['chk_file']
   input_basis_file = j['input_basis_file']
   H_basis_file     = j['H_basis_file']
   At_pp_file       = j['At_pp_file']
   mut_file = j['mutation_file']
   mut_rule = j['mutation_impl']
   cro_rule = j['crossover_impl']

   m_lik = j['mutation_probability'] #Mutation probability is set in 100 in JSON file
   m_mag = j['mutation_magnitude']

   cal_dir          = os.getcwd()
   # Create scratch directories to perform the calculations
   scratch_dir1      = "/lustre/scratch/kruedae/nwchem/1"
   scratch_dir2      = "/lustre/scratch/kruedae/nwchem/2"
   scratch_dir3      = "/lustre/scratch/kruedae/nwchem/3"

   memory     = 32000   # in mb
   xcf        = "pbe0"
#   xhf        = "HFexch" # we will not use HF

   At_title   = "At"
   Atn_title  = "At-"
   Atp_title  = "At+"
   At2p_title  = "At+" #added
   At_m_calc_title = [At2p_title, Atp_title, At_title, Atn_title]
   At2_m_calc_title = ["At2_%f"%(r) for r in R_At2]
   HAt_m_calc_title = ["HAt_%f"%(r) for r in R_HAt]

   At_charge  = 0
   Atn_charge = -1
   Atp_charge = 1
   At2p_charge = 2 #added
   At2_charge = 0
   HAt_charge = 0
   AtI_charge = 0
   At_m_calc_charge = [At2p_charge, Atp_charge, At_charge, Atn_charge]
   At2_m_calc_charge = [At2_charge for l in range(len(R_At2))]
   HAt_m_calc_charge = [HAt_charge for l in range(len(R_HAt))]

   At_mult    = 2
   Atn_mult   = 1
   Atp_mult   = 3
   At2p_mult   = 4 #added
   At2_mult   = 1
   HAt_mult   = 1
   At_m_calc_mult = [At2p_mult, Atp_mult, At_mult, Atn_mult]
   At2_m_calc_mult = [At2_charge for l in range(len(R_At2))]
   HAt_m_calc_mult = [HAt_charge for l in range(len(R_HAt))]

   At_geometry = ["At   0.0   0.0   0.0" for l in range(len(At_m_calc_title))]
   At2_geometry= ["At  0.00000000     0.00000000     0.00000000\n At    0.00000000     0.00000000     %f"%(r) for r in R_At2]
   HAt_geometry= ["H  0.00000000     0.00000000     0.00000000\n At    0.00000000     0.00000000     %f"%(r) for r in R_HAt]
   # basename for the files
   basenameAt  = "nwcat"
   basenameAt2 = "nwcat2"
   basenameHAt = "nwchat"
  
   emin = 1e20
   gens = -1
   inds = -1
   eoriginal = 0.0


   if not os.path.exists(scratch_dir1):
      os.mkdir(scratch_dir1)
   if not os.path.exists(scratch_dir2):
      os.mkdir(scratch_dir2)
   if not os.path.exists(scratch_dir3):
      os.mkdir(scratch_dir3)

   out_file_name = "ga.out"
   fo = open(out_file_name, "w")

   if job_type == 0: # Is job_type a variable to continue interruped calculations or why exists?
      fo.write(" Starting new calculation \n")
      bas_start = basis.read_start_basis(input_basis_file) # I have to change the basis set in the imput file

      fo.write(" -------- Input basis -------- \n")
      for p in range(len(bas_start)):
         v = bas_start[p]
         fo.write("   %s     %18.9f \n"%(v[0], v[1]))
      fo.write(" ------ Input basis end ------ \n\n")

   #
   # read pseudopotential
   #
   fo.write(" Reading pseudopotential from %s \n"%(At_pp_file))
   At_pp = basis.read_pseudopotential(At_pp_file)

   fo.write(" -------- Pseudopotential ----------\n")
   fo.write(At_pp)
   fo.write(" -------- Pseudopotential end ----------\n\n")

   #
   # read H and I basis sets and pseudopotential
   # 
   H_basis = basis.read_atom_basis(H_basis_file)
   #AtI_pp  = basis.read_pseudopotential(AtI_pp_file)  

   #
   # clone input basis into n_ind 
   #
   original_basis_mut = True
   basis_list  = algo.clone(bas_start, n_ind)
   bas_mutated, idxn = algo.mutate(basis_list, original_basis_mut, m_lik, m_mag)

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
      error_At  = np.zeros((n_cur_gen,N_At))
      error_At2  = np.zeros((n_cur_gen,N_At2))
      error_HAt  = np.zeros((n_cur_gen, N_HAt))
      #error_Atp = np.zeros((n_cur_gen))
      #error_Atn = np.zeros((n_cur_gen))

      ips = np.zeros((n_cur_gen))
      eas = np.zeros((n_cur_gen))

      for n, bas in enumerate(bas_mutated):

         # 
         # run NWChem calculations for all individuals in this generation
         # 
         # Atp-Atp-At-At 
         p1 = Process(target=runn.run_many_calc_atom, args=(cal_dir, scratch_dir1, memory, bas, 
                                                  At_pp, At_m_calc_title, At_m_calc_charge, At_m_calc_mult, 
                                                  At_geometry, xhf, basenameAt,))
	 # At2
         p2 = Process(target=runn.run_many_calc_atom, args=(cal_dir, scratch_dir2, memory, bas, 
                                                  At_pp, At2_m_calc_title, At2_m_calc_charge, At2_m_calc_mult, 
                                                  At2_geometry, xhf, basenameAt2,))

	# HAt
         p3 = Process(target=runn.run_many_calc_mol, args=(cal_dir, scratch_dir3, memory, bas, 
                                                  At_pp, HAt_m_calc_title, HAt_m_calc_charge, HAt_m_calc_mult, 
                                                  HAt_geometry, xhf, basenameHAt,H_basis))



         p1.start()
         p2.start()
         p3.start()
         
         p1.join()
         p2.join()
         p3.join()

	 
         E_out_at_m_calc, nwfail1 = runn.get_energy_many_calc(basenameAt, N_At)
         E_out_at2_m_calc, nwfail3 = runn.get_energy_many_calc(basenameAt2, N_At2)
         E_out_hat_m_calc, nwfail2 = runn.get_energy_many_calc(basenameHAt, N_HAt)

	 E_out_at_m_calc = np.array(E_out_at_m_calc)	
	 E_out_at2_m_calc = np.array(E_out_at2_m_calc)	
	 E_out_hat_m_calc = np.array(E_out_hat_m_calc)	
         #--------------------------------------------------------------------------- 

         if nwfail1 or nwfail2 or nwfail3: #or nwfail4 or nwfail5 or nwfail6:
	 	n_failed += 1
	 errors[n] = np.sum(abs(E_out_at_m_calc-E_At_m_calc))
	 errors[n] += np.sum(abs(E_out_at2_m_calc-E_At2))
	 errors[n] += np.sum(abs(E_out_hat_m_calc-E_HAt))
	 
#         error_Atp[n] = Ha_to_eV*np.abs(E_Atp  - E_out_atp)
#         error_Atn[n] = Ha_to_eV*np.abs(E_Atn  - E_out_atn)

	 error_At[n] = Ha_to_eV*abs(E_out_at_m_calc-E_At2)
	 error_At2[n] = Ha_to_eV*abs(E_out_at2_m_calc-E_At2)
	 error_HAt[n] = Ha_to_eV*abs(E_out_hat_m_calc-E_At2)
         #errors[n] += np.abs(E_At2 - E_out_at2)
         #errors[n] += np.abs(E_HAt - E_out_hat)
         #errors[n] += np.abs(E_AtI - E_out_ati)    

         errors[n] *= Ha_to_eV/(N_At+N_At2+N_HAt)

         #ips[n] = Ha_to_eV*(E_out_atp - E_out_at)
         #eas[n] = Ha_to_eV*(E_out_at - E_out_atn)

 

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
      n_promote = int(round(n_rem*1.0/save_n_best)) # Q: Why is a cocient?: Save_n_best must to be a coefficient (10 to 1/10, ...)
      n_tourn   = n_rem - n_promote

      if n_tourn % 2 == 1:
         n_promote += 1
         n_tourn   -= 1     # will make it even

      n_cross = n_ind - n_promote - n_tourn // 2
   
      fo.write(" Directly promoted to the next generation : %d \n"%(n_promote))
      fo.write(" Going into the tournament stage : %d \n"%(n_tourn))

      fo.write(" Best performing basis in this generation gives the error of %10.7f [eV] \n"%(errors[ind[0]]))
      fo.write(" Individual errors: At++ = %10.7f , At+ = %10.7f , At = %10.7f , At- = %10.7f [eV] \n"%(error_At[ind[0],0], error_At[ind[0],1], error_At[ind[0],2], error_At[ind[0],3]))
      fo.write(" Individual errors: At2 = %10.7f , %10.7f , %10.7f , %10.7f , %10.7f [eV] \n"%(error_At2[ind[0],0], error_At2[ind[0],1], error_At2[ind[0],2], error_At2[ind[0],3], error_At2[ind[0],4]))
      fo.write(" Individual errors: HAt = %10.7f , %10.7f , %10.7f , %10.7f , %10.7f [eV] \n"%(error_HAt[ind[0],0], error_HAt[ind[0],1], error_HAt[ind[0],2], error_HAt[ind[0],3], error_HAt[ind[0],4]))
      #fo.write(" IP = %7.5f [eV], EA = %7.5f [eV] \n"%(ips[ind[0]], eas[ind[0]]))
      #fo.write(" Error in IP w.r.t. experiment = %7.5f [eV] \n"%(ips[ind[0]] - At_IP))
      #fo.write(" Error in EA w.r.t. PRA 91, 020501 (2015) = %7.5f [eV] \n"%(eas[ind[0]] - At_EA))
      basis.write_basis(bas_mutated[ind[0]], gen, ind[0]) # Q: Why its writing only the bas_mutated[0] (the best)?
      fo.write(" See file bas_gen%d_id%d.dat \n"%(gen,ind[0]))

       
      err_diff = errors[ind[0]] - errors[ind[n_cur_gen-n_failed-1]]
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
         fo.write(" %d "%(ind_tourn[nx])) #Q: Why is writing the index
      fo.write("\n\n")
      fo.flush()

      #
      # save n best basis sets
      #
      next_gen = [] # Vector with the basis
      for nx in range(n_promote):
         idd = ind[nx] 
         next_gen.append(bas_mutated[idd])#Copy the promoted basis

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
            if rg < fittest_win:# Compares the errors, the lowest wins with probability fittest_win
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
        if cro_rule == 'binary':
           mat_out1, mat_out2 = algo.crossover(next_gen[ind1], next_gen[ind2])
        elif cro_rule == 'float':
           mat_out1, mat_out2 = algo.crossover_float(next_gen[ind1], next_gen[ind2])
        else:
           print (" No crossover rule was given ")
           sys.exit() 
        next_gen.append(mat_out1)
        next_gen.append(mat_out2)
        fo.write(" crossover between %d and %d \n"%(ind1,ind2))
   
      # mutate and go back to the top of the loop
      if mut_rule == "gamma": # m_lik is the probability of mutation, m_mag is a parameter in the mutation value
         m_mag = np.random.gamma(2.0,1.0)
         m_lik = 100.0
         fo.write(" Mutation parameters sampled from Gamma distribution : probability = %7.5f %% magnitude = %7.5f %% \n"%(m_lik, m_mag))
         bas_mutated, idxn = algo.mutate(next_gen, False, m_lik, m_mag)   
      else:
         if os.path.isfile(mut_file):
            m_lik, m_mag = algo.read_mutate(mut_file)
            fo.write(" Mutation parameters read from file %s : probability = %7.5f %% magnitude = %7.5f %% \n"%(mut_file,m_lik, m_mag))
         else:
            fo.write(" Mutation parameters : probability = %7.5f %% magnitude = %7.5f %% \n"%(m_lik, m_mag))  
         bas_mutated, idxn = algo.mutate(next_gen, False, m_lik, m_mag)   
   
      # save checkpoint
      fo.write(" Saving active basis sets to %s \n"%(save_file))
      basis.save(bas_mutated, save_file)

   fo.close()

if __name__ == "__main__":
   main()
