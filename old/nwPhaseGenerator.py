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

W = 7
minN = 7
maxN = 15
Ns = np.arange(minN,maxN+1,1)
M = 0.1
added = False

print("\nGenerating Nanowire Phase (noMagnets = %i --> %i) for M = %1.2f and added: %r..." 
      %(minN,maxN,M,added))

## SOI terms ##
al=.8

print("Also, SO Coupling = %1.1f ..." %(al))

for i in range(np.size(Ns)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=W, noMagnets=Ns[i],
                                  alpha=al, M=M, addedSinu=added
                                  )
    
    ## Phase ##
    pickle.dump(nanowire.phaseTransitions(),
                open("data/phas_" 
                     + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
                     %(W, Ns[i], al, M, int(added))
                     + ".dat", "wb"))
    
print("\nCompleted!")