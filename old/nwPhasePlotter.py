#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np

import pickle
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (8,6)
#plt.rcParams["figure.figsize"] = (15,8)

print("\nPlotting Nanowire Phase transitions...")
W = 7
minN = 17
maxN = 17
Ns = np.arange(minN,maxN+1,1)
M = 0.1
added = False

## SOI terms ##
al=.8

print("Also, SO Coupling = %1.1f ..." %(al))

for i in range(np.size(Ns)):
    ## Phase Data ##
    data = pickle.load(open("data/phas_" 
                            + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
                            %(W, Ns[i], al, M, int(added))
                            + ".dat", "rb"))
    
    
    print("\nPlot for noMagnets = %i" %(Ns[i]))
    
    plt.figure()
    CS = plt.contourf(data["MuSc"], data["B"], data["Phases"], 1, cmap="coolwarm")
    plt.xlabel("Chemical Potential [mu]")
    plt.ylabel("Zeeman Field Strength [B]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Topological State? [Boolean]")
    plt.show()
    
print("\nCompleted!")