import math
import csv
import numpy as np



def gompertz(a,b,c,x):
    '''
    Description: The gompertz function that determines fecundity based on resource level x
    :param a: a is the asymptote of the function 
    :type a: positive floating point number
    :param b: b determines when the function swings up
    :type b: a negative floating point number
    :param c: c determines the rate of approach to the asymptote
    :type c: a negative floating point number
    :param x: the mean resource levels to base reproduction on
    :type x: a single value or list of floating points 
    
    '''
    y = []
    try:
        len(x)
        for s in x:
            y.append(a*math.exp(b*math.exp(c*s)))
        return y
    
    except TypeError:
        return a*math.exp(b*math.exp(c*x))
    
    
def write_ibmdata(LatticeObj):
    '''
    Description: Writes a csv file of all data relevant to our study.  Writes a data file with all the relevant data from the entire simulation.
    right now it writes data on sex ratio, genetic fecundity as calculated by genetics, and then actual fecundity.  Each data item will be written as it's own
    csv file.
    :param LatticeObj: Lattice object from the ibmsimulation library
    :type LatticeObj: Lattice object from the ibmsimulation library
    
    '''
    #extract data from the lattice simulation
    data_toWrite = LatticeObj.output
    open("output.csv","w").close()
    print "Hey I opened the file"
    with open("output.csv","a") as csvfile:
        mywriter = csv.writer(csvfile, delimiter=",")
        for x in data_toWrite:
            mywriter.writerow(x)
    
    
    


