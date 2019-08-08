#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np

import pickle
import matplotlib.pyplot as plt

# L = 20:       0 -> 2
# L = 100:      0 -> 10
# L = 500:      0 -> 50
plt.rcParams["figure.figsize"] = (10,6)

print("\nPlotting Nanowire Data...")
W = 4
L = 70
minPeriod = 0
maxPeriod = .5
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)

for i in range(np.size(periodBs)):
    ## Set up Nanowire Object ##
    print("Plots for periodB = %2.1f" %(periodBs[i]))
    
    ## Spectrum ##
    data = pickle.load(open("zspec" 
                            + "%i.%i.%2.1f" %(W, L, periodBs[i])
                            + ".dat", "rb"))
    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("B")
    plt.ylabel("Energies [t]")
    plt.show()
    
    ## Conductances ##
    data = pickle.load(open("zcond"
                            + "%i.%i.%2.1f" %(W, L, periodBs[i])
                            + ".dat", "rb"))
    plt.figure()
    CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Bias V [t]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.show()

print("\nCompleted!")

# Individual Conductance ##
#data = pickle.load(open("conductance1.dat", "rb"))
#plt.figure()
#plt.plot(data["BiasV"], np.transpose(data["Cond"])[38])
#plt.xlabel("Bias V [t]")
#plt.ylabel("Conductance [e^2/h]")
#plt.show()