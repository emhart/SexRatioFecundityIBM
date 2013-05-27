from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih
import cProfile


#Code profiling...


def test_profile():
    
   
    '''
    Global parameter sets:
    b: Sets the intercept of the gompertz equation, more negative values force the y intercept at 0 longer.
    c: Sets how rapidly the function asymptotes.  more negative values asymptote faster.
    const: sets the threshold reproductive energy needed as a function of total potential fecundity.
    '''
    b = -10
    c = -1
    const = .5

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
    ind_set = {'fr':[1,1],'m_cost':0,'energy':1,'rep_cost': 0 ,'lifespan':1,'fecund_genes':[0,1,2,20],'max_energy':10}



    tmp = L.Lattice(dims = [5,4],Kp = [.01,.02] )

 

    groups = []
    z = []
    for x in range(50):
        indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':1,'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'groupID' : 1,'sex' : np.random.binomial(1,.5,1),'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3])),"max_energy": ind_set["max_energy"]}
        z.append(Ind.individual(**indiv_dict))

    for x in range(20):
        groups.append(G.group([],x,ID=x))

    groups[0] = G.group(z,0,ID=0)

    tmp.groups = groups
    n = 100



    for x in range(n):
        if x%1 == 0:
            print x
    

        tmp.mate(ind_set)
        tmp.disperse(.1,True,24)
        tmp.reproduce()
        tmp.senesce(.05)
        tmp.mutate(0.01)
    
        tmp.regenerate()
        tmp.data_collect()

    ih.write_ibmdata(tmp)
    print "done"

cProfile.run('test_profile()')
