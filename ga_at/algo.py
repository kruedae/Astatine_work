import numpy as np
import copy
import os


def crossover_float(bas1, bas2):
   nf = len(bas1)

   n = 3

   # do a one point crossover here..
   nh = np.random.randint(nf)
   u = np.random.rand()

   bas_out1 = []
   bas_out2 = []

   print (" input bas 1 ",bas1)
   print (" input bas 2 ",bas2)

   print (" changing exponent = ",nh)
  
   for x in range(nf):
      if x == nh:
         di1 = bas1[x] 
         di2 = bas2[x]
         v1 = di1[1]
         v2 = di2[1]
         s1 = di1[0]
         s2 = di2[0] 
         if u < 0.5:
            beta = (2.0*u)**(1.0/(n+1))
         else:
            beta = (1.0/(2.0*(1.0-u)))**(1.0/(n+1))

         real_co1 = 0.5*(v1*(1.0 + beta) + v2*(1.0 - beta))
         real_co2 = 0.5*(v1*(1.0 - beta) + v2*(1.0 + beta))
         tmp1 = []
         tmp1.append(s1)
         tmp1.append(real_co1)
         tmp2 = []
         tmp2.append(s2)
         tmp2.append(real_co2)
         bas_out1.append(tmp1)
         bas_out2.append(tmp2)    
      else:
         bas_out1.append(bas1[x])
         bas_out2.append(bas2[x])

   print (" output bas 1 ",bas_out1)
   print (" output bas 2 ",bas_out2)

   return bas_out1, bas_out2

def read_mutate(mut_file):

   f = open(mut_file,"r")
   f.readline()
   line  = f.readline()
   line2 = line.split()
   m_lik = float(line2[0])
   m_mag = float(line2[1])
   f.close()

   return m_lik, m_mag

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

def mutate(bas_set, mutate_all_but_one, m_lik, m_mag):

   bas_mutated = []
   idxn = -100

   if not mutate_all_but_one:
       for bas_in in bas_set: 
          if np.random.random() < (m_lik/100.0):
             bas_mutated.append(mutate_single(bas_in,m_mag))
          else:
             bas_mutated.append(bas_in)
   else: # mutate_all_but_one makes let a random parameter unmutated
       N = len(bas_set)
       idxn = np.random.randint(N)
       for n in range(N):
          if n == idxn: # if the n-basis is part of the random indxs just make a copy
             bas_mutated.append(bas_set[n])
          else: # if not, the base is mutated with probability m_lik
             if np.random.random() < (m_lik/100.0):
                bas_mutated.append(mutate_single(bas_set[n],m_mag)) # mutate individually the basis set
             else:
                bas_mutated.append(bas_set[n]) 
             
   return bas_mutated, idxn

def mutate_single(bas_in, mag):

   bas_out = copy.deepcopy(bas_in) # use deepcopy to copy the basis set selected to mutate. It is an object?

   N = len(bas_out)
   idx = np.random.randint(N) #a component to mutate is selected
   
   val = bas_out[idx][1]
  
   # allowing to change not larger than 5%
   dv = mag/100.0
   u = np.random.random()
   n = 1
   if u < 0.5:
	delta = (2*u)**(1/(n+1))-1
   if u >= 0.5:
	delta = 1-(2*(1-u))**(1/(n+1))
   #vo = 2.0*dv*val*np.random.random() + val - val*dv # mutation rule performed at the single value
   vo = val+delta*dv # mutation rule performed at the single value

   dd = 100*(vo - val)/val
   #print (" Changing parameter from %10.7f to %10.7f change = %7.5f %%\n"%(val,vo,dd))
   bas_out[idx][1] = vo

   return bas_out
   
