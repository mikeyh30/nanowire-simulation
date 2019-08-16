#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 23:30:48 2019

@author: Domi
"""
import os
import numpy as np
import pickle
import nwObjects
os.system("clear")

# L = 20:       0 -> 2
# L = 100:      0 -> 10
W = 5
L = 100
minPeriod = 0.
maxPeriod = 10.
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)
M = 0.1

print("\nGenerating Nanowire Data (periodB = %2.1f --> %2.1f)..." 
      %(minPeriod,maxPeriod))

for i in range(np.size(periodBs)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=W, length=L, periodB=periodBs[i], M=M)
    
    ## Spectrum ##
    pickle.dump(nanowire.spectrum(bValues=np.linspace(0, 1, 501)),
                open("spec" 
                     + "%i.%i.%2.1f.%1.2f" %(W, L, periodBs[i], M)
                     + ".dat", "wb"))

    ## Conductance ##
    pickle.dump(nanowire.conductances(bValues=np.linspace(0, 1, 101)),
                open("cond"
                     + "%i.%i.%2.1f.%1.2f" %(W, L, periodBs[i], M)
                     + ".dat", "wb"))
    
print("\nCompleted!")