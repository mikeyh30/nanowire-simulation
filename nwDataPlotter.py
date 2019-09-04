#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np

import pickle
import matplotlib.pyplot as plt
#plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Data...")
W = 7
minN = 20
maxN = 20
Ns = np.arange(minN,maxN+1,1)
M = 0.05
added = True

## SOI terms ##
eM=.5
mu=.52
al=.8

for i in range(np.size(Ns)):
    print("\nPlot for noMagnets = %i" %(Ns[i]))
    print("\nSpectrum")
    plt.rcParams["figure.figsize"] = (7,5)
    ## Spectrum ##
    data = pickle.load(open("data/spec_" 
                            + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                            %(W, Ns[i], eM, mu, al, M, int(added))
                            + ".dat", "rb"))
    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Energies [t]")
    plt.show()
    print("Critical value = %1.2f" %(data["CritB"]))
    
    print("\nConductance")
    plt.rcParams["figure.figsize"] = (8,5)
    ## Conductances ##
    data = pickle.load(open("data/cond_" 
                            + "w%i_no%i_eM%1.2f_mu%1.2f_al%1.1f_M%1.2f_added%i" 
                            %(W, Ns[i], eM, mu, al, M, int(added))
                            + ".dat", "rb"))
    plt.figure()
    CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Bias V [t]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.show()
    print("Critical value = %1.2f" %(data["CritB"]))
    
    ## Individual Conductance ##
    plt.rcParams["figure.figsize"] = (7,5)
    index = 30 # 20 & 40
    cond = np.transpose(data["Cond"])
    plt.figure()
    plt.plot(data["BiasV"], cond[index])
    plt.xlabel("Bias V [t]")
    plt.ylabel("Conductance [e^2/h]")
    plt.show()

print("\nCompleted!")