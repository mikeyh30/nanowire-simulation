#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:53:27 2019

@author: Domi
"""

import os
import pickle
import nwObjects
os.system("clear")

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Phase Aid...")
W = 5
N = 5
M = 0.05
added = False

print("\nPlot for noSections = %2.1f" %(N))
nanowire = nwObjects.Nanowire(width=W, noSections=N, M=M, addedSinu=added)
pickle.dump(nanowire.phaseAid(),
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
plt.xlabel("Zeeman field Str, B")
plt.ylabel("Energies [t]")
plt.show()

plt.figure()
plt.plot(data["MuSc"], data["Em"])
plt.xlabel("Chemical Potential, mu")
plt.ylabel("Energies [t]")
plt.show()
    
print("\nCompleted!")