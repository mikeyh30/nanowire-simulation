#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import pickle
import matplotlib.pyplot as plt

n = 8
plt.rcParams["figure.figsize"] = (14,5)

## Spectrum ##
for i in range(n):
    data = pickle.load(open("spectrum" + "%i" %(i) + ".dat", "rb"))

    plt.figure()
    plt.plot(data["B"], data["E"])
    plt.xlabel("B")
    plt.ylabel("Energies [t]")
    plt.show()

## Conductances ##
for i in range(n):
    data = pickle.load(open("conductance" + "%i" %(i) + ".dat", "rb"))

    plt.figure()
    CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
    plt.xlabel("Zeeman Field Strength [B]")
    plt.ylabel("Bias V [t]")
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.show()

## Individuals ##
#data = pickle.load(open("conductance1.dat", "rb"))
#plt.figure()
#plt.plot(data["BiasV"], np.transpose(data["Cond"])[38])
#plt.xlabel("Bias V [t]")
#plt.ylabel("Conductance [e^2/h]")
#plt.show()