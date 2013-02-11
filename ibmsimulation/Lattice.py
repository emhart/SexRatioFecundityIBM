import numpy as np
from collections import Counter
import random as rn
'''
This is my test comment
'''

class Lattice(object):
    '''
    Testing docs here just the doc of my class
    masterData holds all the data for a simulation
    '''
    def __init__(self,dims,rate,init_energy,max_energy,groups=[],size=0):
        self.rate = rate
        self.energy = init_energy
        self.occupied = [0]*(dims[0]*dims[1])
        self.xmax = dims[0]
        self.max_energy = max_energy
        self.groups = groups
        self.size = (dims[0]*dims[1])
        self.timeStep = 0
        #This needs to be the same size right now this has columns for time step, colony id, group size, actualFecund mean, var, genetic fecund mean, var 
        self.output = np.array([1,2,3,4,5,6,7])
        
        
        
    def distance(self,pos):
        '''
        Oh hey here's more documentanion
        param: pos - a vector position
        '''
        coord = divmod(pos[0],self.xmax)
        if coord[1] > 0:
            p1 = np.array([coord[1],coord[0]+1])
        else:
            p1=np.array([coord[0],self.xmax])
        coord = divmod(pos[1],self.xmax)
        if coord[1] > 0:
            p2 = np.array([coord[1],coord[0]+1])
        else:
            p2=np.array([coord[0],self.xmax])

        return np.linalg.norm(p2-p1)
    
    def regenerate(self):
        '''
        Description: Regenerate patch energy by looping through each patch and adding energy
        Right now I will leave it as a normal variate that can be set.  This allows for environmental
        stochasticity.  
        :param mu: The mean amount of resource that should be added at each time step back to the patch
        :type mu: positive float
        :param sigma: The variance in the normal draw made, corresponds to the total amount of environmental noise
        :type sigma: 
        '''
        self.energy = self.energy+self.rate
        self.energy[self.energy > self.max_energy] = self.max_energy
    
    def disperse(self,d_p):
        '''
        Description: Randomly sort through each group, colonize new patches randomly:
        '''
        ## Set up index for randomized dispersal
        disp_index = rn.sample(range(self.size),self.size)
        for x in disp_index:
            to_disp = self.groups[x].disperse()
            for i in to_disp:
                #### Tune dispersal with a simple binomial for now....
                if np.random.binomial(1,d_p) == 1: 
                    self.groups[np.random.randint(0,self.size)].indivs.append(i)

    def mutate(self):
        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].mutate() 
                
                
    def reproduce(self):
        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].reproduce()

                
    def mate(self,b,c,const,ind_set):
        off_num = []
        for i in range(self.size):
            if self.groups[i]:
                off_num.append(self.groups[i].mate(b,c,const,ind_set))
        #remove zero values 
        off_num = filter (lambda a: a != 0, off_num)
 
        
    def forage(self):
        for i in range(self.size):
            if self.groups[i]:
                self.energy = self.groups[i].forage(self.energy)
                
    def senesce(self,mort_prob):

        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].senesce(mort_prob)    
                    
    def data_collect(self):
        '''
        Description:  This will collect all the relevant data from different colonies and store it in 
        :modifies: self.geneticFecund
        :modifies: self.actcualFecund
        :modifies: self.groupSize
        '''

        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].report()
                mean_genF = self.groups[i].geneticFecund
                var_genF = self.groups[i].genetic_FVar
                mean_actF = self.groups[i].actFecund
                var_actF = self.groups[i].actual_FVar
                #meanForage = np.append(meanForage,self.groups[i].meanForage)
                #meanEn = np.append(meanEn,self.groups[i].meanEnergy)
                self.groups[i].resize()
                pop_sizes = self.groups[i].size
                curID = self.groups[i].ID
                self.output = np.vstack((self.output,[self.timeStep,curID,pop_sizes,mean_genF,var_genF,mean_actF,var_actF]))
            if not self.groups[i]:
                self.output = np.vstack((self.output[self.timeStep,self.groups[i].ID,0,0,0,0,0]))
        
        self.timeStep += 1                
                
        
    def mean_size(self):
        
        for i in self.groups:
            i.resize()
            pop_sizes = np.append(pop_sizes,i.size) 
        return pop_sizes

        
        