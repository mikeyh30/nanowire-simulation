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
N = 10
M = 0.1
added = False

print("\nPlot for noMagnets = %2.1f" %(N))
nanowire = nwObjects.Nanowire(width=W, noMagnets=N, M=M, addedSinu=added)
pickle.dump(nanowire.phaseAid(bValues=np.linspace(.0, .18, 37),
                              muValues=np.linspace(.08, .2, 49)
                              ),
            open("aide" 
                 + "w%i.no%i.m%1.2f.added%i" %(W, N, M, int(added))
                 + ".dat", "wb"))

## Phase Data ##
data = pickle.load(open("aide" 
                        + "w%i.no%i.m%1.2f.added%i" %(W, N, M, int(added))
                        + ".dat", "rb"))
## Plot individual ##
plt.figure()
plt.plot(data["B"], data["Eb"])
plt.xlabel("Zeeman field Str [B]")
plt.ylabel("Energies [t]")
plt.ylim([-.12, .12])
plt.show()

plt.figure()
plt.plot(data["MuSc"], data["Em"])
plt.xlabel("Chemical Potential [mu]")
plt.ylabel("Energies [t]")
plt.ylim([-.12, .12])
plt.show()
    
print("\nCompleted!")