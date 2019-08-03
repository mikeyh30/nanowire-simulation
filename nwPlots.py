#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 01:47:15 2019

@author: Domi
"""
import os
os.system("clear")

import numpy as np
import pickle
import nwObjects
import matplotlib.pyplot as plt

## Set up Nanowire Objects (Restrict number of object in a single run) ##
nanowire0 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=0)
nanowire1 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=0.5)
nanowire2 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=1.0)
nanowire3 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=1.5)
nanowire4 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=2.0)
nanowire5 = nwObjects.Nanowire(width=5, length=50, barrierLen=1, periodB=2.0)
nanowire6 = nwObjects.Nanowire(width=5, length=50, barrierLen=1, periodB=10.5)
nanowire7 = nwObjects.Nanowire(width=5, length=100, barrierLen=1, periodB=10.5)

plt.rcParams["figure.figsize"] = (14,5)

## Spectrum ##
print("Calculating Spectrum(s)...\n")
data = nanowire0.spectrum()
pickle.dump(data, open("spectrum0.dat", "wb"))
plt.figure()
plt.plot(data["B"], data["E"])
plt.xlabel("B")
plt.ylabel("Energies")
plt.show()

data = nanowire1.spectrum()
pickle.dump(data, open("spectrum1.dat", "wb"))

data = nanowire2.spectrum()
pickle.dump(data, open("spectrum2.dat", "wb"))

pickle.dump(nanowire3.spectrum(), open("spectrum3.dat", "wb"))
pickle.dump(nanowire4.spectrum(), open("spectrum4.dat", "wb"))
pickle.dump(nanowire5.spectrum(), open("spectrum5.dat", "wb"))
pickle.dump(nanowire6.spectrum(), open("spectrum6.dat", "wb"))
pickle.dump(nanowire7.spectrum(), open("spectrum7.dat", "wb"))
print("Completed!\n")

## Conductances ##
print("Calculating Conductance(s)...\n")
data = nanowire0.conductances(bValues=np.linspace(0, 1.0, 101))
pickle.dump(data, open("conductance0.dat", "wb"))

plt.figure()
CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
plt.xlabel("Zeeman Field Strength [B]")
plt.ylabel("Bias V [t]")
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('Conductance')
plt.show()

data = nanowire1.conductances(bValues=np.linspace(0, 1.0, 101))
pickle.dump(data, open("conductance1.dat", "wb"))

data = nanowire2.conductances(bValues=np.linspace(0, 1.0, 101))
pickle.dump(data, open("conductance2.dat", "wb"))

pickle.dump(nanowire3.conductances(bValues=np.linspace(0, 1.0, 101)),
            open("conductance3.dat", "wb"))
pickle.dump(nanowire4.conductances(bValues=np.linspace(0, 1.0, 101)),
            open("conductance4.dat", "wb"))
pickle.dump(nanowire5.conductances(bValues=np.linspace(0, 1.0, 101)),
            open("conductance5.dat", "wb"))
pickle.dump(nanowire6.conductances(bValues=np.linspace(0, 1.0, 101)),
            open("conductance6.dat", "wb"))
pickle.dump(nanowire7.conductances(bValues=np.linspace(0, 1.0, 101)),
            open("conductance7.dat", "wb"))
print("Completed!\n")