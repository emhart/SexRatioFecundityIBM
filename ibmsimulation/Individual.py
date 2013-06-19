import numpy as np
from collections import Counter
from ibmsimulation import ibm_help as ih


'''
@author: Edmund Hart

'''


class individual(object):
    '''
    Description: A class that holds all the information for a given individual.
    :param age: the current age
    :type age: a scalar integer
    :param sex: the sex of the individual 1 is male, 0 is female
    :type sex: scalar integer
    '''
    def __init__(self,forage_rate,m_cost,rep_cost,lifespan,age=0,energy=[],sex=0,fecund_genes=[],groupID = 0,starv_time = 0, mutation_r = .05,pregnant=[],max_energy=[]):
        self.age = age
        self.energy = energy
        self.sex = sex
        self.fecund_genes = fecund_genes
        self.forage_rate = forage_rate
        self.m_cost = m_cost
        self.rep_cost = rep_cost
        self.lifespan = lifespan
        self.groupID = groupID
        self.starv_time= starv_time
        self.pregnant = pregnant
        self.food_hist = [0,0]
        self.max_energy = max_energy
        
    def meiosis(self,myChromo,sizes):
        '''
        Description: A function that returns a single chromosome from meisosis. Allows for crossing over using 
        :param myChromo: The chromosome to have undergo meisosis
        :type myChromo: a 2 x L numpy array where L is the number of loci in the chromosome. 
        :returns: a 1 x L numpy array that will be used in gamete creation
        '''
            
        c_size = len(myChromo[0])

        my_randint = np.random.randint(sizes[0][c_size-1],sizes[1][c_size-1],1)
        as_bin = np.array([int(d) for d in bin(my_randint)[2:]],dtype=bool)
        bin_comp = np.array([int(d) for d in bin(np.uint32(~ my_randint))[(len(bin(np.uint32(~ my_randint)))-len(as_bin)):34]],dtype=bool)


        myChromo = np.append(myChromo,np.array([np.append(myChromo[0][as_bin],myChromo[1][bin_comp]),np.append(myChromo[1][as_bin],myChromo[0][bin_comp])]),axis=0)
        
        return myChromo[np.random.randint(0,4,1)]
    
    
    def forage(self,pos,energy):
        '''
        Description:  This allows individuals to forage.  It transfers energy to the organism if 
        it is available and also extracts a cost.  An organism may not necessarily get any energy, but always incurs 
        at maintenance cost. 
        
        :param pos: The position of the colony
        :type pos: scalar int
        :param c_lattice: The underlying lattice of the model
        :type c_lattice: an object of class Lattice
        :returns: int -- the new energy of the patch after the individual has foraged
        :modifies: self.starv_time
        
        '''
        pot_energy = energy[pos]
        if pot_energy > self.forage_rate:
            self.energy = self.forage_rate + self.energy 
            return pot_energy - self.forage_rate
        
        elif pot_energy < self.forage_rate and pot_energy > 0: 
            self.energy = self.energy + pot_energy 
            return 0

        elif pot_energy <= 0:
            self.energy = self.energy 
            return 0
                      
    def create_gametes(self,sizes):
        fecundity_gamete = self.meiosis(self.fecund_genes,sizes)
        return fecundity_gamete
    def mutate(self,rate):
        '''
        Description: Mutates an individuals genome, set up for a continuous genome
        :param rate: The probability that an individual mutates at each time step (this should be very low)
        :type rate: Floating point between 0 and 1
        
        '''
        
        if np.random.binomial(1,rate,1) == 1:
            #select a chromosome
            chr = np.random.choice([0,1],1)
            loci = np.random.choice(range(len(self.fecund_genes[chr,])),1)
            self.fecund_genes[chr,loci] = np.random.uniform(0,1,1)

            
                
            
            
        
