import numpy as np
import copy
import os

def crossover(bas1, bas2):
   nf = len(bas1)

   # do a one point crossover here..
   nh = np.random.randint(1,nf-1)

   bas_out1 = []
   bas_out2 = []
   for nx in range(nh):
      bas_out1.append(bas1[nx])
      bas_out2.append(bas2[nx])
   for nx in range(nh,nf):
      bas_out1.append(bas2[nx])
      bas_out2.append(bas1[nx])

   return bas_out1, bas_out2     

def clone(bas_in, nclone):
   
   bas_set_out = []
   for n in range(nclone):
      bas_set_out.append(bas_in)

   return bas_set_out

def mutate(bas_set, mutate_all_but_one):

   m_lik = 100.0
   m_mag = 1.0
 
   # read mutation parameters
   if os.path.isfile('mutate.txt'):
      f = open("mutate.txt","r")
      f.readline()
      line  = f.readline()
      line2 = line.split()
      m_lik = float(line2[0])
      m_mag = float(line2[1])
      f.close()
 
   bas_mutated = []
   idxn = -100
 
   if not mutate_all_but_one:
       for bas_in in bas_set: 
          if np.random.random() < (m_lik/100.0):
             bas_mutated.append(mutate_single(bas_in,m_mag))
          else:
             bas_mutated.append(bas_in)
   else:
       N = len(bas_set)
       idxn = np.random.randint(N)
       for n in range(N):
          if n == idxn:
             bas_mutated.append(bas_set[n])
          else:
             if np.random.random() < (m_lik/100.0):
                bas_mutated.append(mutate_single(bas_set[n],m_mag))
             else:
                bas_mutated.append(bas_set[n]) 
             
   return bas_mutated, idxn

def mutate_single(bas_in, mag):

   bas_out = copy.deepcopy(bas_in)

   N = len(bas_out)
   idx = np.random.randint(N)   
   
   val = bas_out[idx][1]
  
   # allowing to change not larger than 5%
   dv = mag/100.0

   vo = 2.0*dv*val*np.random.random() + val - val*dv

   dd = 100*(vo - val)/val
   #print (" Changing parameter from %10.7f to %10.7f change = %7.5f %%\n"%(val,vo,dd))
   bas_out[idx][1] = vo

   return bas_out
   
