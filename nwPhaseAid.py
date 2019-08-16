#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:53:27 2019

@author: Domi
"""

import os
import numpy as np
import pickle
import nwObjects
os.system("clear")

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Phase Aid...")
# L = 20:       0 -> 2
# L = 100:      0 -> 10
W = 5
L = 100
periodB = 0.5
M = 0.05

print("\nPlot for periodB = %2.1f" %(periodB))

nanowire = nwObjects.Nanowire(width=W, length=L, periodB=periodB, M=M)
pickle.dump(nanowire.phaseAid(),
            open("aide" 
                 + "%i.%i.%2.1f.%1.2f" %(W, L, periodB, M)
                 + ".dat", "wb"))

## Phase Data ##
data = pickle.load(open("aide" 
                        + "%i.%i.%2.1f.%1.2f" %(W, L, periodB, M)
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