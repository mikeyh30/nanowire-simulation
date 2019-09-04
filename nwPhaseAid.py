#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:53:27 2019

@author: Domi
"""

import os
import numpy as np
import nwObjects
os.system("clear")

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (5,4)

print("\nPlotting Nanowire Phase Aid...")
W = 7
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
data = nanowire.phaseAid(bValues=np.linspace(.0, .2, 41),
                         muValues=np.linspace(.1, .8, 71) # add 24
                         )

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