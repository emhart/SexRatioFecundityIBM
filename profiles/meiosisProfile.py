    
import numpy as np
import cProfile

s = ((2**0,2**1,2**2,2**3,2**4,2**5,2**6,2**7,2**8,2**9,2**10,2**11,2**12,2**13,2**14,2**15,2**16,2**17),((2**1-1,2**2-1,2**3-1,2**4-1,2**5-1,2**6-1,2**7-1,2**8-1,2**9-1,2**10-1,2**11-1,2**12-1,2**13-1,2**14-1,2**15-1,2**16-1,2**17-1,2**18-1)))

def meiosis(myChromo,sizes):

    c_size = len(myChromo[0])

    my_randint = np.random.randint(sizes[0][c_size-1],sizes[1][c_size-1],1)
    #as_bin = np.array([int(d) for d in bin(my_randint)[2:]],dtype=bool)
    #bin_comp = np.array([int(d) for d in bin(np.uint32(~ my_randint))[(len(bin(np.uint32(~ my_randint)))-len(as_bin)):34]],dtype=bool)
    as_bin = np.array(map(int,bin(my_randint)[2:]),dtype=bool)
    bin_comp = np.array(map(int,bin(np.uint32(~ my_randint))[(len(bin(np.uint32(~ my_randint)))-len(as_bin)):34]),dtype=bool)

    tmp = np.append(myChromo,np.array([np.append(myChromo[0][as_bin],myChromo[1][bin_comp]),np.append(myChromo[1][as_bin],myChromo[0][bin_comp])]),axis=0)
        

def m_test():
    for x in range(100000):
        x = np.random.uniform(0,1,[2,5])
        meiosis(x,s)

cProfile.run('m_test()')