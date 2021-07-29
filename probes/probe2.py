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


#------------ Reference data ----------
Ha_to_eV = 27.211
E_At  = -263.0812543    # SO HF/uncontracted-cc-pV5Z-PP atomic calc. starting from anion orbitals
E_Atp = -262.7419035    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
E_At2p = -262.0950157    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
E_Atn = -263.1663307    # SO HF/uncontracted-cc-pV5Z-PP anion calc. starting from the default guess
# Put the values of distance and energy
R_At2 = [2.97842834, 2.77842834, 2.87842834, 3.07842834, 3.17842834 ]
E_At2 = [-526.1925101, -526.1867749, -526.191208, -526.1916841, -526.1894469]
R_HAt = [1.5237857, 1.6237857, 1.7237857, 1.8237857, 1.9237857]
E_HAt = [-263.6547155, -263.6666387, -263.670009, -263.6677744, -263.6619657]
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
# CREATE A NEW DIRECTORY FOR EACH CALCULATION (14)
scratch_dir1      = "/lustre/scratch/kruedae/nwchem/1"
scratch_dir2      = "/lustre/scratch/kruedae/nwchem/2"
scratch_dir3      = "/lustre/scratch/akanane/nwchem/3"
scratch_dir4      = "/lustre/scratch/akanane/nwchem/4"
scratch_dir5      = "/lustre/scratch/akanane/nwchem/5"
scratch_dir6      = "/lustre/scratch/akanane/nwchem/6"
scratch_dir7      = "/lustre/scratch/akanane/nwchem/7"
scratch_dir8      = "/lustre/scratch/akanane/nwchem/8"
scratch_dir9      = "/lustre/scratch/akanane/nwchem/9"
scratch_dir10      = "/lustre/scratch/akanane/nwchem/10"
scratch_dir11      = "/lustre/scratch/akanane/nwchem/11"

memory     = 32000   # in mb
xcf        = "pbe0"
#   xhf        = "HFexch" # we will not use HF

At_title   = "At"
Atn_title  = "At-"
Atp_title  = "At+"
At2p_title  = "At+" #added
At2_title  = "At2"
HAt_title  = "HAt"

At_charge  = 0
Atn_charge = -1
Atp_charge = 1
At2p_charge = 2 #added
At2_charge = 0
HAt_charge = 0
AtI_charge = 0

At_mult    = 2
Atn_mult   = 1
Atp_mult   = 3
At2p_mult   = 4 #added
At2_mult   = 1
HAt_mult   = 1
At_geometry = "At   0.0   0.0   0.0"
At2_geometry= ["At  0.00000000     0.00000000     0.00000000\n At    0.00000000     0.00000000     %f"%(r) for r in R_At2]
HAt_geometry= ["H  0.00000000     0.00000000     0.00000000\n At    0.00000000     0.00000000     %f"%(r) for r in R_HAt]
# basename for the files
basenameAt  = "nwc"
basenameAtp = "nwcp"
basenameAt2p = "nwc2p"
basenameAtn = "nwcn"
basenameAt2 = ["nwcd%s"%(r) for r in np.arange(len(R_At2))]
basenameHAt = ["nwchat%s"%(r) for r in np.arange(len(R_HAt))]

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
if not os.path.exists(scratch_dir4):
os.mkdir(scratch_dir4)
if not os.path.exists(scratch_dir5):
os.mkdir(scratch_dir5)
if not os.path.exists(scratch_dir6):
os.mkdir(scratch_dir6)
if not os.path.exists(scratch_dir7):
os.mkdir(scratch_dir7)
if not os.path.exists(scratch_dir8):
os.mkdir(scratch_dir8)
if not os.path.exists(scratch_dir9):
os.mkdir(scratch_dir9)
if not os.path.exists(scratch_dir10):
os.mkdir(scratch_dir10)
if not os.path.exists(scratch_dir11):
os.mkdir(scratch_dir11)

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
error_At2  = np.zeros((n_cur_gen, R_At2))
error_HAt  = np.zeros((n_cur_gen, R_HAt))
#error_Atp = np.zeros((n_cur_gen))
#error_Atn = np.zeros((n_cur_gen))

ips = np.zeros((n_cur_gen))
eas = np.zeros((n_cur_gen))

for n, bas in enumerate(bas_mutated):

 # 
 # run NWChem calculations for all individuals in this generation
 # 
 # At2 
 p1 = Process(target=runn.run_atom, args=(cal_dir, scratch_dir1, memory, bas, 
					  At_pp, At2_title, At2_charge, A2_mult, 
					  At2_geometry[0], xhf, basenameAt2[0],))
