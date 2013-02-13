'''
Created on Jul 19, 2012


@author: emh
A working file to test out my python classes, will be junked when the full simulation is run.

'''
from ibmsimulation import Lattice as L
from ibmsimulation import Individual as Ind
from ibmsimulation import Group as G
import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih

#Code profiling...



ID_stack = range(1000000)
ID_stack.reverse()

b = -10
c = -10
const = 2

#Set the parameters for all individual
ind_set = {'fr':[.1,.8],'m_cost':.3,'energy':20,'rep_cost':.1,'lifespan':20,'rep_thresh': 20,'fecund_genes':[0,1,2,4]}



tmp = L.Lattice(dims = [3,3],rate = np.random.uniform(100,120,25),init_energy = [100]*25, max_energy = 1000 )



groups = []
z = []
for x in range(3):
    indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':ind_set["energy"],'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'ID' : ID_stack.pop(),'groupID' : 1,'sex' : np.random.binomial(1,.5,1),'rep_thresh': ind_set["rep_thresh"],'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3]))}
    z.append(Ind.individual(**indiv_dict))

for x in range(9):
    groups.append(G.group([],x,ID=x))

groups[0] = G.group(z,8,ID=0)

tmp.groups = groups
n = 40



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


'''
n = 2000
size_dat = np.array([.1]*n)
aveF_dat = np.zeros(n)
for x in range(n):
    if x%50 == 0:
        print x
    groups[10].forage(tmp)
    groups[10].reproduce(indiv_dict)
    groups[10].resize()
    groups[10].sr_count()
    tmp.regenerate()
    groups[10].mutate()
    groups[10].senesce()
    aveF_dat[x] = groups[10].meanFecund
    size_dat[x] = groups[10].size 
#print groups[10].indivs[4].lifespan 
#print groups[10].indivs[4].age


groups[10].age_count()
#
plt.plot(size_dat)
print groups[10].indivs[5].starv_time
plt.show()
r_dat= np.log(size_dat[1:(n-1)]) - np.log(size_dat[0:(n-2)])
plt.scatter(size_dat[1:(n-2)],r_dat[1:(n-2)])
plt.show()
print groups[10].indivs[30].ID
plt.plot(aveF_dat)
plt.ylabel("Mean fecundity per individual within a group")
plt.xlabel("Time step")
plt.title("Evolution of fecundity")
plt.show()


groups[10].sr_count()
print groups[10].sex_r


#groups[10].sr_count()

##it = Ind.individual(forage_rate=2,m_cost=1,energy=5,rep_cost=20,rep_thresh=40,fecund_genes=np.random.uniform(0,1,(2,16)))

'''

