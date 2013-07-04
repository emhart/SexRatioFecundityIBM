import numpy as np
from collections import Counter
from ibmsimulation import ibm_help as ih
import re


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
    def __init__(self,sex=0,fecund_genes=[],groupID = 0, mutation_r = .05,pregnant=[], age = 0,lifespan = 1):
        self.sex = sex
        self.fecund_genes = fecund_genes
        self.groupID = groupID
        self.pregnant = pregnant
        self.age = age
        self.lifespan = lifespan

        
    def meiosis(self,myChromo):
        '''
        Description: Crossing over happens at only one spot now.  
        :param myChromo: The chromosome to have undergo meisosis
        :type myChromo: a 2 x L numpy array where L is the number of loci in the chromosome. 
        :returns: a 1 x L numpy array that will be used in gamete creation
        '''
        #choose a random cut point
        cp = np.random.randint(1,15,1)
        
        #Convert to binary
        c1 = '{0:0>32b}'.format(myChromo[0])
        c2 = '{0:0>32b}'.format(myChromo[1])
        #cross over
        cx1 = c1[cp:] + c2[:cp]
        cx2 = c2[cp:] + c1[:cp]
        #return to integers
        return np.random.choice([int(c1,2),int(c2,2),int(cx1,2),int(cx2,2)])
    
    

    def create_gametes(self):
        fecundity_gamete = self.meiosis(self.fecund_genes)
        return fecundity_gamete
    
    def mutate(self,rate):
        '''
        Description: Mutates an individuals genome, set up for a continuous genome
        :param rate: The probability that an individual mutates at each time step (this should be very low)
        :type rate: Floating point between 0 and 1
        
        '''
        
        if np.random.binomial(1,rate,1) == 1:
            ## Create binary string from gene
            chr_ind = np.random.choice([0,1],1)           
            chr_str = '{0:0>32b}'.format(int(self.fecund_genes[chr_ind]))
            #get location of 1 and 0
            one = re.finditer("1",chr_str)
            one_ind = [n.start() for n in one]
            zero = re.finditer("0",chr_str)
            zero_ind = [n.start() for n in zero]            
            chr_str = list(chr_str)
            ###Now check for mutation possibilities
            if zero_ind and one_ind:
                chr_str[np.random.choice(one_ind,1)] = "0"
                chr_str[np.random.choice(zero_ind,1)] = "1"
                chr_str = ''.join(chr_str)
                self.fecund_genes[chr_ind] = int(chr_str,2)
            elif zero_ind and not one_ind:
                chr_str[np.random.choice(zero_ind,1)] = "1"
                chr_str = ''.join(chr_str)
                self.fecund_genes[chr_ind] = int(chr_str,2)                
            elif one_ind and not zero_ind:
                chr_str[np.random.choice(one_ind,1)] = "0"
                chr_str = ''.join(chr_str)
                self.fecund_genes[chr_ind] = int(chr_str,2)                
            
                
            
            
        
