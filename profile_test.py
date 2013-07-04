from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih
import cProfile


#Code profiling...


def test_profile():
    
 




    ind_set = {'fecund_genes':[1,2**16]}



    tmp = L.Lattice(dims = [5,4],Kp = [.001,.02] )

    n_patch = 20
    groups = []

    for i in range(n_patch):
        z = []
        for x in range(20):
            indiv_dict = {'groupID' : n_patch,'sex' : np.random.binomial(1,.5,1),'fecund_genes':np.random.randint(1,2**3,2)}
            z.append(Ind.individual(**indiv_dict))
    
        groups.append(G.group(z,i,ID=i))

    tmp.groups = groups
    n = 100

    for x in range(n):
        if x%1 == 0:
            print x
    
    
        tmp.mate(ind_set)
    #tmp.disperse(.1,True,24)
        tmp.colonize(disp_size = 10,  migration_p  = 0 ,migration_frac=.1)
 
        tmp.reproduce()

        tmp.senesce(.05)

        tmp.mutate(0.1)
        tmp.regenerate()
    
        tmp.data_collect()

    ih.write_ibmdata(tmp)
    print "done"



cProfile.run('test_profile()')
