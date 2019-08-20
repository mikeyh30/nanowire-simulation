# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:54:11 2019

@author: dk1713
"""

import os
import numpy as np
import pickle
import nwObjects
os.system("clear")

W = 5
minN = 3
maxN = 17
Ns = np.arange(minN,maxN+1,2)
M = 0.05
added = True

print("\nGenerating Nanowire Phase (noSections = %i --> %i) and added: %r..." 
      %(minN,maxN, added))

for i in range(np.size(Ns)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=W, noSections=Ns[i], M=M, addedSinu=added)
    
    ## Phase ##
    pickle.dump(nanowire.phaseTransition(),
                open("data/phas" 
                     + "w%i.no%i.m%1.2f.added%i" %(W, Ns[i], M, int(added))
                     + ".dat", "wb"))
    
print("\nCompleted!")