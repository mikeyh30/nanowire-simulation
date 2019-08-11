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
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (15,7)
os.system("clear")
# L = 20:       0 -> 2
# L = 100:      0 -> 10
# L = 500:      0 -> 50
W = 5
L = 100
minPeriod = 0
maxPeriod = 2.
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)

print("\nGenerating Nanowire Data (periodB = %2.1f --> %2.1f)..." 
      %(minPeriod,maxPeriod))
for i in range(np.size(periodBs)):
    ## Set up Nanowire Object ##
    nanowire = nwObjects.Nanowire(width=W, length=L, periodB=periodBs[i])
    
    ## Spectrum ##
    pickle.dump(nanowire.spectrum(bValues=np.linspace(0, 1, 51)),
                open("tspec" 
                     + "%i.%i.%2.1f" %(W, L, periodBs[i])
                     + ".dat", "wb"))
    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("B")
    plt.ylabel("Energies [t]")
    plt.show()

    ## Conductances ##
    pickle.dump(nanowire.conductances(bValues=np.linspace(0, 1, 51)),
                open("tcond"
                     + "%i.%i.%2.1f" %(W, L, periodBs[i])
                     + ".dat", "wb"))
print("\nCompleted!")