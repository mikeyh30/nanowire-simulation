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
plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Data...")
W = 5
L = 100
minPeriod = 0
maxPeriod = .5
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)

critB = []

for i in range(np.size(periodBs)):
    ## Set up Nanowire Object ##
    print("Plots for periodB = %2.1f" %(periodBs[i]))
    
    ## Spectrum ##
    data = pickle.load(open("tspec" 
                            + "%i.%i.%2.1f" %(W, L, periodBs[i])
                            + ".dat", "rb"))
    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("B")
    plt.ylabel("Energies [t]")
    plt.show()
    
    ## Conductances ##
    data = pickle.load(open("tcond"
                            + "%i.%i.%2.1f" %(W, L, periodBs[i])
                            + ".dat", "rb"))
    plt.figure()
    CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Bias V [t]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.show()
    
    print("Critical value = %1.2f" %(data["CritB"]))
    critB.append(1/data["CritB"])

print("\nCompleted!")

critB = np.multiply(critB,1/critB[0])
#d = np.multiply(np.add(periodBs,0.5),2)
#d[0] = 0

plt.figure()
plt.plot(periodBs, critB)
plt.xlabel("d")
plt.ylabel("Critical B in Uniform/Critical B in Sinusoidal")
plt.show()

# Individual Conductance ##
#data = pickle.load(open("conductance1.dat", "rb"))
#plt.figure()
#plt.plot(data["BiasV"], np.transpose(data["Cond"])[38])
#plt.xlabel("Bias V [t]")
#plt.ylabel("Conductance [e^2/h]")
#plt.show()