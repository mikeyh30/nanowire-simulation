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
# L = 20:       0 -> 2
# L = 100:      0 -> 10
W = 5
L = 100
minPeriod = 0
maxPeriod = 10.
periodBs = np.arange(minPeriod,maxPeriod+.5,.5)
M = 0.1

print("\nPlotting Nanowire Phases..." )

ratio = []
critical = []

for i in range(np.size(periodBs)):
#    print("\nPlot for periodB = %2.1f" %(periodBs[i]))
    
    ## Phase Data ##
    data = pickle.load(open("phas" 
                            + "%i.%i.%2.1f.%1.2f" %(W, L, periodBs[i], M)
                            + ".dat", "rb"))
    ## Plot individual ##
#    plt.figure()
#    plt.plot(data["MuSc"], data["CritB"])
#    plt.xlabel("Chemical Potential, mu")
#    plt.ylabel("Critical Point, B")
#    plt.show()
    
    ## Ratio: for muSc = (indices), sinusoidal/uniform ##
    ratio.append(data["CritB"][40])
    critical.append(data["CritB"])
    
print("\nCompleted!")

## Plotting Ratio for muSc = (indices) ##
ratio = np.multiply(ratio,1/ratio[0])
plt.figure()
plt.plot(periodBs, ratio, 'b-')
plt.xlabel("d")
plt.ylabel("Critical Point with Sinusoidal fields/Critical Point with Uniform fields")
plt.show()

## 3D Phase Diagram ##
fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(periodBs, data["MuSc"])
critical = np.transpose(critical)

# Plot the surface.
surf = ax.plot_surface(X, Y, critical, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
cbar = fig.colorbar(surf, shrink=0.5, aspect=5)
cbar.ax.set_ylabel("Critical Point, B")
ax.set_xlabel("d")
ax.set_ylabel("Chemical Potential, mu")
ax.set_zlabel("Critical Point, B")
plt.show()