'''
Created on Jul 19, 2012


@author: emh
A working file to test out my python classes, will be junked when the full simulation is run.

'''

### Imports

from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih



'''
Set the parameters for all individuals
fr: The upper and lower bounds of the feeding rate parameter, drawn from a uniform distribution
energy: The starting energy of a newly born individual
rep_cost: The energetic cost of reproduction per individual offspring    
lifespan: The number of time steps an organism can live.
rep_thresh: The energetic threshold that an organism needs to reach.
fecund_genes: A list of four numbers. Positions [0,1] are the upper and lower bounds of a uniform distribution
and [2] is the number of chromosomes, usually 2, and position [3] is the length of each chromosome.

'''
ind_set = {'fecund_genes':[1,2**16]}



tmp = L.Lattice(dims = [10,5],Kp = [-4.7,.7  ] )

n_patch = 50
groups = []

for i in range(n_patch):
    z = []
    for x in range(20):
        indiv_dict = {'groupID' : n_patch,'sex' : np.random.binomial(1,.5,1),'fecund_genes':np.random.randint(1,2**2,2)}
        z.append(Ind.individual(**indiv_dict))
    
    groups.append(G.group(z,i,ID=i))

tmp.groups = groups
n = 1000

for x in range(n):
    if x%1 == 0:
        print x
    
    
    tmp.mate(ind_set)
    tmp.colonize(disp_size = 10,  migration_p  = 0 ,migration_frac=.2,mult=1)
 
    tmp.reproduce()

    tmp.senesce(.05)

    tmp.mutate(0.01)
    tmp.regenerate()
    
    tmp.data_collect()

ih.write_ibmdata(tmp)
print "done"



