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

plt.rcParams["figure.figsize"] = (4.5,3)

print("\nPlotting Nanowire Phase transitions...")
W = 7
minN = 7
maxN = 20
Ns = np.arange(minN,maxN+1,1)
M = 0.1

## SOI terms ##
al=.0

ratio0 = []
ratio1 = []
index = 0

critical0 = []
critical1 = []


added = False
for i in range(np.size(Ns)):
    ## Phase Data ##
    data = pickle.load(open("data/crit_" 
                            + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
                            %(W, Ns[i], al, M, int(added))
                            + ".dat", "rb"))
    ## Plot individual ##
    if Ns[i] == 10:
#        plt.rcParams["figure.figsize"] = (4,4)
        print("\nPlot for noMagnets = %i" %(Ns[i]))
        plt.figure()
#        plt.plot(data["MuSc"][0:30], data["CritB"][0:30])
        plt.plot(data["MuSc"], data["CritB"])
        plt.xlabel("Chemical Potential [mu]")
        plt.ylabel("Critical Point [B]")
#        plt.ylim([0,.3])
        plt.show()
    
    ## Ratio: for muSc = (indices), sinusoidal/uniform ##
    ratio0.append(data["CritB"][index])
    critical0.append(data["CritB"])
    
print("\nCompleted!")

added = True
for i in range(np.size(Ns)):
    ## Phase Data ##
    data = pickle.load(open("data/crit_" 
                            + "w%i_no%i_al%1.1f_M%1.2f_added%i" 
                            %(W, Ns[i], al, M, int(added))
                            + ".dat", "rb"))
    ## Plot individual ##
    if Ns[i] == 10:
        print("\nPlot for noMagnets = %i" %(Ns[i]))
        plt.figure()
        plt.plot(data["MuSc"], data["CritB"])
        plt.xlabel("Chemical Potential [mu]")
        plt.ylabel("Critical Point [B]")
        plt.show()
    
    ## Ratio: for muSc = (indices), sinusoidal/uniform ##
    ratio1.append(data["CritB"][index])
    critical1.append(data["CritB"])
    
print("\nCompleted!")

## Plotting Ratio for muSc = (index) ##
plt.figure()
plt.plot(Ns, np.array(ratio1)/ratio0, 'b-')
plt.xlabel("No of Magnets")
plt.ylabel("Critical Point Ratio")
plt.show()

plt.rcParams["figure.figsize"] = (9,6)
## 3D Phase Diagram ##
print("Uniform Field")
fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(Ns, data["MuSc"])
critical0 = np.transpose(critical0)

# Plot the surface.
surf = ax.plot_surface(X, Y, critical0, cmap=cm.coolwarm, rcount=100,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
cbar.ax.set_ylabel("Critical Point [B]")
ax.set_xlabel("No of Magnets")
ax.set_ylabel("Chemical Potential [mu]")
ax.set_zlabel("Critical Point [B]")
ax.set_zlim(0, .3)
plt.show()