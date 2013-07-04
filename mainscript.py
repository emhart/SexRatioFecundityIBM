#!/usr/bin/python
### Imports

from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import numpy as np
from ibmsimulation import ibm_help as ih
import sys
import re



def read_args():
    args_file = sys.argv[1]
    d = {}
    with open(args_file) as f:
        for line in f:
            (key, val) = line.split()
            d[key] = float(val)
    return(d)


def main():
    '''
        Description:  Runs a full simulation of our spider model.  Requires a command line input of a file with parameters.
        Below is a full description of parameters.  They must be named as described, separated by spaces.
        ngens: The number of generations in the simulation
        x_dim: The x dimensions of the lattice that we are simulating over.
        y_dim: The y dimension of the lattice we are simulating over.
        min_pop: The minimum population size expressed as c in 1/c = pop size, e.g. min_pop = 0.01 means a minimum population size of 100
        max_pop: The maximum population size expressed as c in 1/c = pop size, e.g. min_pop = 0.01 means a minimum population size of 100
    '''

    #extract parameters
    param_dict = read_args()

    #Set global parameters
    x_dim = int(param_dict['x_dim'])
    y_dim = int(param_dict['y_dim'])
    min_pop = param_dict['min_pop']
    max_pop = param_dict['max_pop']
    ngens = int(param_dict['ngens'])
    min_col_size = int(param_dict['min_col_size'])


    ind_set = {'fr':[1,1],'m_cost':0,'energy':1,'rep_cost': 0 ,'lifespan':1,'fecund_genes':[0,1,2,20],'max_energy':10}



    tmp = L.Lattice(dims = [x_dim,y_dim],Kp = [max_pop,min_pop] )

    n_patch = x_dim * y_dim

    groups = []

    for i in range(n_patch):
        z = []
        for x in range(20):
            indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':1,'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'groupID' : 1,'sex' : np.random.binomial(1,.5,1),'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3])),"max_energy": ind_set["max_energy"]}
            z.append(Ind.individual(**indiv_dict))
    
        groups.append(G.group(z,i,ID=i))

    tmp.groups = groups
    n = ngens



    for x in range(n):
        if x%1 == 0:
            print x
    

        tmp.mate(ind_set)
        #tmp.disperse(.1,True,24)
        tmp.colonize(min_col_size,.5,.001)
        tmp.reproduce()
        tmp.senesce(.05)
        tmp.mutate(0.001)
        #tmp.forage()
        tmp.regenerate()
        tmp.data_collect()

    fName = "data/output_min_pop_"+str(min_pop)+"_min_col_size_"+str(min_col_size)+".csv"
    ih.write_ibmdata(tmp,fName)
    print "done"

if __name__ == "__main__":
    main()

