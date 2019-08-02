#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 01:47:15 2019

@author: Domi
"""
import os
os.system("clear")

import pickle
import nwObjects
import matplotlib.pyplot as plt

nanowire0 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=0)
nanowire1 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=0.5)
nanowire2 = nwObjects.Nanowire(width=5, length=20, barrierLen=1, periodB=1)
plt.rcParams["figure.figsize"] = (8,5)

## Spectrum ##
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

## Conductances ##
data = nanowire0.conductances()
pickle.dump(data, open("conductance0.dat", "wb"))

plt.figure()
CS = plt.contourf(data["B"], data["BiasV"], data["Cond"], 100, cmap="viridis")
plt.xlabel("Zeeman Field Strength [B]")
plt.ylabel("Bias V [t]")
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('Conductance')
plt.show()

data = nanowire1.conductances()
pickle.dump(data, open("conductance1.dat", "wb"))

data = nanowire2.conductances()
pickle.dump(data, open("conductance2.dat", "wb"))