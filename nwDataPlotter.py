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

print("\nPlotting Nanowire Data...")
W = 7
minN = 7
maxN = 20
Ns = np.arange(minN,maxN+1,1)
M = 0.1
added = True

for i in range(np.size(Ns)):
    print("\nPlot for noMagnets = %i" %(Ns[i]))
    print("\nSpectrum")
    plt.rcParams["figure.figsize"] = (4.5,3)
    ## Spectrum ##
    data = pickle.load(open("data/spec" 
                            + "w%i.no%i.m%1.2f.added%i" %(W, Ns[i], M, int(added))
                            + ".dat", "rb"))
    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Energies [t]")
    plt.show()
    print("Critical value = %1.2f" %(data["CritB"]))
    
    print("\nConductance")
    plt.rcParams["figure.figsize"] = (5,3)
    ## Conductances ##
    data = pickle.load(open("data/cond" 
                            + "w%i.no%i.m%1.2f.added%i" %(W, Ns[i], M, int(added))
                            + ".dat", "rb"))
    plt.figure()
    CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Bias V [t]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.show()
    print("Critical value = %1.2f" %(data["CritB"]))

print("\nCompleted!")