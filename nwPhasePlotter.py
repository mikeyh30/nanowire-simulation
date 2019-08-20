#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 04:02:18 2019

@author: Domi
"""

import numpy as np

import pickle
from mpl_toolkits.mplot3d import Axes3D # not an error
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams["figure.figsize"] = (15,7)

print("\nPlotting Nanowire Phase transitions...")
W = 5
minN = 3
maxN = 15
Ns = np.arange(minN,maxN+1,2)
M = 0.05

ratio0 = []
ratio1 = []
index = 0

critical0 = []
critical1 = []

added = False
for i in range(np.size(Ns)):
    print("\nPlot for noSections = %i" %(Ns[i]))
    
    ## Phase Data ##
    data = pickle.load(open("data/phas" 
                            + "w%i.no%i.m%1.2f.added%i" %(W, Ns[i], M, int(added))
                            + ".dat", "rb"))
    ## Plot individual ##
    plt.figure()
    plt.plot(data["MuSc"], data["CritB0"])
    plt.xlabel("Chemical Potential, mu")
    plt.ylabel("Critical Point, B")
    plt.show()
    
    ## Ratio: for muSc = (indices), sinusoidal/uniform ##
    ratio0.append(data["CritB0"][index])
    critical0.append(data["CritB0"])
    
print("\nCompleted!")

added = True
for i in range(np.size(Ns)):
    print("\nPlot for noSections = %i" %(Ns[i]))
    
    ## Phase Data ##
    data = pickle.load(open("data/phas" 
                            + "w%i.no%i.m%1.2f.added%i" %(W, Ns[i], M, int(added))
                            + ".dat", "rb"))
    ## Plot individual ##
    plt.figure()
    plt.plot(data["MuSc"], data["CritB0"])
    plt.xlabel("Chemical Potential, mu")
    plt.ylabel("Critical Point, B")
    plt.show()
    
    ## Ratio: for muSc = (indices), sinusoidal/uniform ##
    ratio1.append(data["CritB0"][index])
    critical1.append(data["CritB0"])
    
print("\nCompleted!")

## Plotting Ratio for muSc = (index) ##
plt.figure()
plt.plot(Ns, ratio1/ratio0, 'b-')
plt.xlabel("d")
plt.ylabel("Critical Point with Sinusoidal fields/Critical Point with Uniform fields")
plt.show()

## 3D Phase Diagram ##
print("Uniform Field")
fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(Ns, data["MuSc"])
critical0 = np.transpose(critical0)

# Plot the surface.
surf = ax.plot_surface(X, Y, critical0, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
cbar.ax.set_ylabel("Critical Point, B")
ax.set_xlabel("d")
ax.set_ylabel("Chemical Potential, mu")
ax.set_zlabel("Critical Point, B")
plt.show()

## 3D Phase Diagram ##
fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(Ns, data["MuSc"])
critical1 = np.transpose(critical1)

# Plot the surface.
print("Sinu Field")
surf = ax.plot_surface(X, Y, critical1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
cbar.ax.set_ylabel("Critical Point, B")
ax.set_xlabel("d")
ax.set_ylabel("Chemical Potential, mu")
ax.set_zlabel("Critical Point, B")
plt.show()