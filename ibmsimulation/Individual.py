import numpy as np
from collections import Counter
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
    def __init__(self,forage_rate,m_cost,rep_cost,rep_thresh,lifespan,age=0,energy=20,sex=0,fecund_genes=[],ID=0,parentID_m=[],parentID_f=[],groupID = 0,starv_time = 0, mutation_r = .05,pregnant=[]):
        self.age = age
        self.energy = energy
        self.sex = sex
        self.fecund_genes = fecund_genes
        self.forage_rate = forage_rate
        self.m_cost = m_cost
        self.rep_cost = rep_cost
        self.rep_thresh = rep_thresh
        self.lifespan = lifespan
        self.ID = ID
        self.parentID_m = parentID_m
        self.parentID_f = parentID_f
        self.groupID = groupID
        self.starv_time= starv_time
        self.mutation_r = mutation_r
        self.pregnant = pregnant
        self.food_hist = [0,0,0,0,0]
        
    def meiosis(self,myChromo):
        '''
        Description: A function that returns a single chromosome from meisosis.  Allows for crossing over at up to two points with a probability that you set.  Crossing ever events happen
        with probability bin_cx.  Needs to be run on each chromosome that you want to undergo meisosis.
        :param myChromo: The chromosome to have undergo meisosis
        :type myChromo: a 2 x L numpy array where L is the number of loci in the chromosome. 
        :returns: a 1 x L numpy array that will be used in gamete creation
        '''
        bin_cx = .5
        outer_chromosomes = np.copy(myChromo)

        inner_chromosomes = np.copy(myChromo)
        ic_copy = np.copy(inner_chromosomes)
        # calculate the probability of a crossing over event
        cx_prob = np.random.binomial(2,bin_cx)
        if cx_prob == 1:
            cp = np.random.randint(0,len(outer_chromosomes[1]))
                      

            inner_chromosomes[0,range(cp)] = ic_copy[1,range(cp)]
            inner_chromosomes[1,range(cp)] = ic_copy[0,range(cp)]
            to_ret = np.random.binomial(1,.5,size=2)
            cells = [outer_chromosomes,inner_chromosomes]
            return cells[to_ret[0]][to_ret[1]]
        if cx_prob == 2:
            cp = np.random.randint(0,len(outer_chromosomes[1]),size=2)
            cp = sorted(cp)
            
            inner_chromosomes[0,range(cp[0])] = ic_copy[1,range(cp[0])]
            inner_chromosomes[1,range(cp[0])] = ic_copy[0,range(cp[0])]
            inner_chromosomes[0,range(cp[1],len(outer_chromosomes[1]))] = ic_copy[1,range(cp[1],len(outer_chromosomes[1]))]
            inner_chromosomes[1,range(cp[1],len(outer_chromosomes[1]))] = ic_copy[0,range(cp[1],len(outer_chromosomes[1]))]

            to_ret = np.random.binomial(1,.5,size=2)
            cells = [outer_chromosomes,inner_chromosomes]
            return cells[to_ret[0]][to_ret[1]]
        else:
            return outer_chromosomes[np.random.binomial(1,.5)]
    
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
            self.energy = self.forage_rate + self.energy - self.m_cost
            self.food_hist.pop()
            self.food_hist.insert(0,self.forage_rate)
            
            self.starv_time = 0
            return pot_energy - self.forage_rate
        elif pot_energy < self.forage_rate and pot_energy > 0: 
            self.energy = self.energy + pot_energy - self.m_cost
            self.food_hist.pop()
            self.food_hist.insert(0,pot_energy)
            
            self.starv_time +=1
            return 0
        elif pot_energy <= 0:
            self.energy = self.energy - self.m_cost
            self.food_hist.pop()
            self.food_hist.insert(0,0)
            
            self.starv_time +=1
            return 0
                      
    def create_gametes(self):
        fecundity_gamete = self.meiosis(self.fecund_genes)
        return fecundity_gamete
    
    def mutate(self):
        '''
        Description: Mutates an individuals genome, set up for a continuous genome
        
        '''
        
        if np.random.binomial(1,self.mutation_r,1) == 1:
            #select a chromosome
            chr = np.random.choice([0,1],1)
            loci = np.random.choice(range(len(self.fecund_genes[chr,])),1)
            self.fecund_genes[chr,loci] = self.fecund_genes[chr,loci] + np.random.choice([-1,1],1)*np.random.uniform(0,1,1)
            if self.fecund_genes[chr,loci] > 1:
                self.fecund_genes[chr,loci] = 1
            if self.fecund_genes[chr,loci] < 0:
                self.fecund_genes[chr,loci] = 0
            
                
            
            
        
