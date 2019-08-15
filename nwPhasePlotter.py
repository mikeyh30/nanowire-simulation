#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np

import pickle
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Phase transitions...")
# L = 20:       0 -> 2
# L = 100:      0 -> 10
W = 5
L = 100
minPeriod = 0
maxPeriod = 0.
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)
M = 0.1

print("\nGenerating Nanowire Data (periodB = %2.1f --> %2.1f)..." 
      %(minPeriod,maxPeriod))

for i in range(np.size(periodBs)):
    print("Plots for periodB = %2.1f" %(periodBs[i]))
    
    ## Phase ##
    data = pickle.load(open("phas" 
                            + "%i.%i.%2.1f.%1.2f" %(W, L, periodBs[i], M)
                            + ".dat", "rb"))
    plt.figure()
    plt.plot(data["MuSc"], data["CritB"])
    plt.xlabel("Chemical Potential, mu")
    plt.ylabel("Critical Point, B")
    plt.show()
    
print("\nCompleted!")