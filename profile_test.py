from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih
import cProfile


#Code profiling...


def test_profile():
    
 



    ind_set = {'fr':[1,1],'m_cost':0,'energy':1,'rep_cost': 0 ,'lifespan':1,'fecund_genes':[0,.5,2,25],'max_energy':10}



    tmp = L.Lattice(dims = [5,2],Kp = [.01,.02] )

    n_patch = 10
    groups = []

    for i in range(n_patch):
        z = []
        for x in range(20):
            indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':1,'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'groupID' : 1,'sex' : np.random.binomial(1,.5,1),'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3])),"max_energy": ind_set["max_energy"]}
            z.append(Ind.individual(**indiv_dict))
    
        groups.append(G.group(z,i,ID=i))

    tmp.groups = groups
    n = 500



    for x in range(n):
        if x%1 == 0:
            print x
    

        tmp.mate(ind_set)
        #tmp.disperse(.1,True,24)
        tmp.colonize(10000,1,.01)
        tmp.reproduce()
        tmp.senesce(.05)
        tmp.mutate(0.001)
        #tmp.forage()
        tmp.regenerate()
        tmp.data_collect()

    ih.write_ibmdata(tmp)
    print "done"

cProfile.run('test_profile()')
