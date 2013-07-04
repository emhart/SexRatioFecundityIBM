import numpy as np
from collections import Counter
import random as rn
from itertools import chain
import math
'''
This is my test comment
'''

class Lattice(object):
    '''
    Testing docs here just the doc of my class
    masterData holds all the data for a simulation
    '''
    def __init__(self,dims,Kp,groups=[],size=0):
        self.occupied = np.array([0]*(dims[0]*dims[1])) 
        self.xmax = dims[0]
        self.groups = groups
        self.size = (dims[0]*dims[1])
        self.ID = range(self.size)
        self.groupID = range(self.size)
        self.timeStep = 0
        #This needs to be the same size right now this has columns for time step, colony id, group size, actualFecund mean, var, genetic fecund mean, var 
        self.output = np.array([1,2,3,4,5,6,7,8,9])
        # Just a t
        self.sizesTup = ((2**0,2**1,2**2,2**3,2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16,2**17,2**18,2**19,2**20,2**21,2**22,2**23,2**24,2**25),((2**1-1,2**2-1,2**3-1,2**4-1,2**5-1,2**6-1,2**7-1,2**8-1,2**9-1,2**10-1,2**11-1,2**12-1,2**13-1,2**14-1,2**15-1,2**16-1,2**17-1,2**18-1,2**19-1,2**20-1,2**21-1,2**22-1,2**23-1,2**24-1,2**25-1)))
        self.Kp = Kp       
        self.K = np.random.lognormal(Kp[0],Kp[1],self.size)
        self.Km1 = self.K
        
        
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
        Description: Carrying capacity of a patch is defined as 1/K, and should satisfy the condition  0 < x < 1
        regenerating involves drawing a random uniform number and resetting K.
        :param Kp: a list giving the lower and upper bounds of a uniform distribution
        :type Kp: list of size 2 of floating point numbers between 0 and 1
        :modifies K: Adjusts the carrying capacity of the patch

        '''
        self.Km1 = self.K
        self.K = np.random.lognormal(self.Kp[0],self.Kp[1],self.size)
    
    
                    

    def dispersal(self,mult):
        '''
        Description: Returns a list of potential dispersers
        :returns dispersers: a dictionary where patch ID is the Key, and the values are a list of dispersers 
        
        '''
        dispersers = {}
        for ind, g in enumerate(self.groups):
            dispersers[self.ID[ind]] = g.disperse(self.Km1[ind],mult)
        return dispersers
        
    def colonize(self,disp_size = 2, migration_p = 0, migration_frac = .001,mult = 1):
        '''
         Description: a colonization and migration function.  It performs two tasks.  First it allows the colonization of empty sites with a given propagule size.  
        Also it allows for migration via a migration probability
        '''
        self.set_occupied()
        to_disp = self.dispersal(mult)
        ## Get the proportion for each dispersal
        disp_numbers = np.array(map(len,to_disp.values()),dtype = np.float64)
        disp_probs = disp_numbers / sum(disp_numbers)
        # get the propagule number
        prop_number = disp_numbers // disp_size 
        
        empty_sites = np.where(self.occupied == 0)[0]  
        
        if empty_sites.size > 0 and sum(prop_number) > 0:
            ### Draw colonizers
            col_number = np.random.multinomial(empty_sites.size,disp_probs)
            
            #Truncate to fit the number of possible colonizers
            col_number[prop_number < col_number] = prop_number[prop_number < col_number]
            for ID, propagules in enumerate(col_number):
                for colonizer in range(propagules):
                    # Take a random sample from the right spot in the dictionary
                    prop = rn.sample(to_disp[ID],disp_size)
                    #remove the disperser
                    for i in prop:
                        to_disp[ID].remove(i)
                    ### now append them to their new location ###
                    
                    # draw a random site
                    new_site = rn.choice(empty_sites)
                    self.groups[new_site].indivs.extend(prop)
                    self.groupID[new_site] = ID
                    
                    #remove that site from the list
                    self.occupied[new_site] = 1  
                    empty_sites = np.where(self.occupied == 0)[0]  
                disp_index = rn.sample(range(self.size),self.size)
        
        ### migration component
        disp_index = rn.sample(range(self.size),self.size)
        ## flatten our dispersers
        migrate_pool =  list(chain.from_iterable(to_disp.values()))
        for x in disp_index:
            #### Tune dispersal with a simple binomial for now....
            if np.random.binomial(1,migration_p) == 1: 
                num_migrate = int(math.floor(self.groups[x].size * migration_frac))
                #now pop it off 
                if len(migrate_pool) > num_migrate:
                    migrators = rn.sample(migrate_pool,num_migrate)
                    self.groups[x].indivs.extend(migrators)
                    
                    


    def mutate(self,rate):
        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].mutate(rate) 
                
                
    def reproduce(self):
        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].reproduce()

                
    def mate(self,ind_set):

        for i in range(self.size):
            if self.groups[i]:
                self.groups[i].mate(self.K[i],ind_set)
       
        
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
                p_ID = self.ID[i]
                g_ID = self.groupID[i]

                self.output = np.vstack((self.output,[self.timeStep,curID,p_ID,g_ID,pop_sizes,mean_genF,var_genF,mean_actF,var_actF]))
            if not self.groups[i]:
                self.output = np.vstack((self.output[self.timeStep,self.groups[i].ID,0,0,0,0,0,0,0]))
        
        self.timeStep += 1                
                
        
    def mean_size(self):
        
        for i in self.groups:
            i.resize()
            pop_sizes = np.append(pop_sizes,i.size) 
        return pop_sizes
    
    def set_occupied(self):
        '''
        Description: Checks each patch to see if its occupied or not.  Sets the appropriate status, occupied or not.
        :modifies occupied: sets the occupied flag to 1 or 0 depending on migration or extinction
        
        '''
        #occ = []
        for i in range(self.size):
            if self.groups[i].size > 0:
                self.occupied[i] = 1
                #self.groupID[i] = self.groups[i].ID
                
            else:
                self.occupied[i] = 0
                self.groupID[i] = 9999
                
            
            
        
        
        