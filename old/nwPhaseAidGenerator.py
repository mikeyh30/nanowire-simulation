#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:53:27 2019

@author: Domi
"""

import os
import pickle
import numpy as np
import nwObjects
os.system("clear")

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (4.5,3)

print("\nPlotting Nanowire Phase Aid...")
W = 7
minN = 12
maxN = 15
N = 20
M = 0.1
added = False

## SOI terms ##
al=.4

print("Also, SO Coupling = %1.1f ..." %(al))

## Set up Nanowire Object ##
nanowire = nwObjects.Nanowire(width=W, noMagnets=N,
                              alpha=al, M=M, addedSinu=added
                              )

## Phase ##
pickle.dump(nanowire.phaseAid(bValues=np.linspace(.0, .15, 31),
                              muValues=np.linspace(.1, .34, 49) # add 24
                              ),
            open("data/aide_" 
                 + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
                 %(W, N, al, M, int(added))
                 + ".dat", "wb"))
    
print("\nCompleted!")

#for i in [True, False]: # added?
#    for j in [0.05, 0.1]: # M
#        for k in [0.0, 0.4, 0.8] # alpha
#            ## Set up Nanowire Object ##
#            nanowire = nwObjects.Nanowire(width=W, noMagnets=N,
#                                          alpha=k, M=j, addedSinu=i
#                                          )
#            
#            ## Phase ##
#            pickle.dump(nanowire.phaseAid(bValues=np.linspace(.0, .18, 37),
#                                          muValues=np.linspace(.08, .2, 49)
#                                          ),
#                        open("data/aide_" 
#                             + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
#                             %(W, Ns[i], al, M, int(added))
#                             + ".dat", "wb"))
#    
#print("\nCompleted!")