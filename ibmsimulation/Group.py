from collections import Counter
import random as ran
import numpy as np
from ibmsimulation import Individual as Ind
from ibmsimulation import ibm_help 

class group(object):
    '''
    Description: A class group which defines all actions that happen on the group level.  Each square in the lattice is
    a group.  
    
    
    '''
    def __init__(self,indivs,pos,ages={},sex_ratio={},meanFecund = 0,ID = 0,ID_stack=range(1000000),**indiv_dict):
        self.indivs=indivs
        self.size = len(indivs)
        self.ages = ages
        self.pos = pos   
        self.ID = ID       
        self.IDmax = 0    
        self.ID_stack = ID_stack
        ###Reporter data structures. 
        #Actual fecundity
        self.actFecund = []
        ## Genetic Fecundity
        self.geneticFecund = []
        #variables for variance
        self.genetic_FVar = []
        self.actual_FVar = []
        #sex ratio
        self.sex_r = []
        ##individual associated with each trait
        #Individual ID
        #self.IDs
        #Paternal ID
        #self.pmIDs
        #Maternal ID
        #self.pfIDs
        
        
    def age_count(self):
        '''
        Description:  Count the number of individuals in a group and update the size attribute
                      Also counts the sex ratio at the same time so no need to call that method separately.
        '''
        
        ###extract data from colony object into a vector
        ### and return a dictionary unique identifiers
        agec = []
        for i in range(self.size):
            agec.append(self.indivs[i].age)
        agecount = Counter(agec)
        self.ages = dict(agecount)
        
    ##### generic reporter function that will set group level data based on individual level properties.      
    def report(self):
    ###extract data from colony object into a vector
        '''
        Description:  This reports data of interest from individuals for each time step.  A helper function 
        write_ibmdata() will extract information from group objects.
        
        :modifies: self.geneticFecund
        :modifies: self.genetic_FVar
        :modifies: self.sex_r
        '''
        sr = []
        tf = []
        #fr = []
        #enlev = []
        
        
        for i in range(self.size):
            sr.append(self.indivs[i].sex[0])
            tf.append(self.indivs[i].fecund_genes.sum())
           # fr.append(self.indivs[i].food_hist[0])
           # enlev.append(self.indivs[i].energy)
        
        #Set values
        srcount = Counter(sr)
        tf = np.array(tf)
        
        ### Set reporter objects        
        #self.sex_r = append(srcount.values())
        #fr = np.array(fr)
        #enlev = np.array(enlev)
        #self.meanForage = fr.mean()
        #self.meanEnergy = enlev.mean()

        self.geneticFecund = tf.mean()
        self.genetic_FVar = tf.var()
        
        
        
        
    def resize(self):
        '''
        Description:  Count the number of individuals in a group and update the size attribute
                      Also counts the sex ratio at the same time so no need to call that method separately.
        '''
        self.size = len(self.indivs)
        
    def forage(self,energy):
        '''
        Have every individual forage on the patch
        '''    
        
        for x in ran.sample(range(self.size),self.size):
            energy[self.pos] = self.indivs[x].forage(self.pos,energy)
            self.indivs[x].age = self.indivs[x].age + 1
        return energy
            
    def mate(self,b,c,const,ind_set):
        ''',
         Description:  Sexual reproduction within the group.  Random mating of all females.  Males can mate
        multiple times.
        
         
        '''
        repro_count_f = np.array([]) #create an index of the number of females that reproduce
        repro_count_m = np.array([]) #create an index of the number of males that reproduce
        tmp_IDstack = self.ID_stack
        
        off_num = np.array([]) #The actual number of off spring reproduced, returned 
        
        for x in range(self.size):
            if self.indivs[x].sex == 1 and self.indivs[x].energy > (np.sum(self.indivs[x].fecund_genes) / const):
                repro_count_f = np.append(repro_count_f,x)
            #if self.indivs[x].sex == 0 and self.indivs[x].energy > self.indivs[x].rep_thresh:
                #repro_count_m = np.append(repro_count_m,x)
            
            # In this version all males can reproduce no matter what
            if self.indivs[x].sex == 0:
                repro_count_m = np.append(repro_count_m,x)


        
        repro_count_f = repro_count_f.astype(np.int64)
        repro_count_m = repro_count_m.astype(np.int64)
        if repro_count_m.sum() > 0:
            for x in repro_count_f:
                babies = []
                mating_m = np.random.choice(repro_count_m,1)
                ### Set the max fecundity
                max_fecund = (np.sum(self.indivs[x].fecund_genes) + np.sum(self.indivs[mating_m].fecund_genes)) /2
                ### lamda for the poison distribution
                lam = ibm_help.gompertz(max_fecund,b,c,sum(self.indivs[x].food_hist) / len(self.indivs[x].food_hist))
                
                offspring_size = np.random.poisson(lam)
                self.indivs[x].energy = self.indivs[x].energy - (offspring_size*self.indivs[x].rep_cost)
                
                # Keep track of the average number of offspring when females reproduce
                off_num = np.append(off_num,offspring_size)
               
                for i in range(offspring_size):
                    ID = tmp_IDstack.pop()
                    indiv_dict = {'forage_rate':np.random.uniform(ind_set["fr"][0],ind_set["fr"][1]),'m_cost':ind_set["m_cost"],'energy':ind_set["energy"],'rep_cost':ind_set["rep_cost"],'lifespan':ind_set["lifespan"],'ID' : ID,'groupID' : self.ID,'sex' : np.random.binomial(1,.5,1),'rep_thresh': ind_set["rep_thresh"],'fecund_genes':np.random.uniform(ind_set["fecund_genes"][0],ind_set["fecund_genes"][1],(ind_set["fecund_genes"][2],ind_set["fecund_genes"][3]))}

                    indiv_dict['fecund_genes'] = np.array([self.indivs[x].create_gametes(),self.indivs[mating_m].create_gametes()])
                
                    indiv_dict['parentID_f'] = self.indivs[x].ID
                    indiv_dict['parentID_m'] = self.indivs[mating_m].ID
                    babies.append(indiv_dict)
            
                self.indivs[x].pregnant = babies
        
        self.resize()   
        self.ID_stack = tmp_IDstack
        if off_num.any():
            self.actFecund = off_num.mean()
            self.actual_FVar = off_num.var()
        elif not off_num.any():
            self.actFecund = 0
            self.actual_FVar = 0
        
    def reproduce(self):
        
        self.resize()
        
        for x in range(self.size):
            if self.indivs[x].pregnant:
                for i in self.indivs[x].pregnant:
                    self.indivs.append(Ind.individual(**i))
                self.indivs[x].pregnant = []
         
   
    def senesce(self,mort_prob):
        '''
        Description: Checks if any individual has zero energy or less, and if they've reached their life span.
        If so they are removed from the population.  
        
        :param mort_prob: the background mortality probability.  Causes random mortality in the group.
        :type mort_prob: floating point, between 0 and 1    
        :modifies: self.indivs
        
        '''
        i_copy = list(self.indivs)
        for x in i_copy: #use this len calculationto ensure the proper size is used
            
            if np.random.binomial(1,mort_prob) == 1:
                self.indivs.remove(x)
                
            
            elif x.age > x.lifespan and x.energy > 0:
                self.indivs.remove(x)
                #print "I died of old age"
            
            elif x.energy <= 0:
                self.indivs.remove(x)
                #print "I died of no food"
            
            
        ### resize when done
        self.resize()    

    def mutate(self):
        for x in self.indivs:
            x.mutate()
                   
    def disperse(self):
        '''
        Description: checks the starving time on each individual and removes those individuals from the current group.  
        :modifies: self.indivs
        :returns: a list of individuals that are dispersing
        '''
        dispersers = []
        i_copy = list(self.indivs)
        for x in i_copy:
            if x.starv_time >= 2:
                dispersers.append(x)
                self.indivs.remove(x)
        return dispersers


