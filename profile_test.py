from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih
import cProfile


#Code profiling...


def test_profile():
    
    ID_stack = range(1000000)
    ID_stack.reverse()

    b = -1
    c = -10
    const = 2

#Set the parameters for all individual
    ind_set = {'fr':[.1,.8],'m_cost':.3,'energy':20,'rep_cost':.1,'lifespan':20,'rep_thresh': 20,'fecund_genes':[0,1,2,10]}



    tmp = L.Lattice(dims = [3,3],rate = np.random.uniform(100,120,25),init_energy = [100]*25, max_energy = 1000 )



    groups = []
    z = []
    for x in range(30):
        indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':ind_set["energy"],'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'ID' : ID_stack.pop(),'groupID' : 1,'sex' : np.random.binomial(1,.5,1),'rep_thresh': ind_set["rep_thresh"],'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3]))}
        z.append(Ind.individual(**indiv_dict))

    for x in range(9):
        groups.append(G.group([],x,ID=x))


    groups[0] = G.group(z,8,ID=0)

    tmp.groups = groups
    n = 100



    for x in range(n):
        if x%20 == 0:
            print x
        tmp.forage()
        tmp.regenerate()
        tmp.mutate()
        tmp.mate(b,c,const,ind_set)
        tmp.disperse(1)
        tmp.reproduce()
        tmp.senesce(0.05)
        tmp.data_collect()

    ih.write_ibmdata(tmp)
    print "done"
    

cProfile.run('test_profile()')
