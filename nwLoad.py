#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import pickle
import matplotlib.pyplot as plt

## Spectrum ##

data = pickle.load(open("spectrum0.dat", "rb"))

plt.figure()
plt.plot(data["B"], data["E"])
plt.xlabel("B")
plt.ylabel("Energies")
plt.show()

data = pickle.load(open("spectrum1.dat", "rb"))

plt.figure()
plt.plot(data["B"], data["E"])
plt.xlabel("B")
plt.ylabel("Energies")
plt.show()

data = pickle.load(open("spectrum2.dat", "rb"))

plt.figure()
plt.plot(data["B"], data["E"])
plt.xlabel("B")
plt.ylabel("Energies")
plt.show()

## Conductances ##

data = pickle.load(open("conductance0.dat", "rb"))

plt.figure()
CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
plt.xlabel("Zeeman Field Strength [B]")
plt.ylabel("Bias V [t]")
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('Conductance')
plt.show()

data = pickle.load(open("conductance1.dat", "rb"))

plt.figure()
CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
plt.xlabel("Zeeman Field Strength [B]")
plt.ylabel("Bias V [t]")
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('Conductance')
plt.show()

data = pickle.load(open("conductance2.dat", "rb"))

plt.figure()
CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
plt.xlabel("Zeeman Field Strength [B]")
plt.ylabel("Bias V [t]")
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('Conductance')
plt.show()

## Individuals
#data = pickle.load(open("conductance1.dat", "rb"))
#plt.figure()
#plt.plot(data["BiasV"], np.transpose(data["Cond"])[38])
#plt.xlabel("Bias V [t]")
#plt.ylabel("Conductance [e^2/h]")
#plt.show()