#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 03:16:59 2019

@author: Domi
"""

import kwant
import tinyarray as ta
#import numpy as np
import matplotlib.pyplot as plt

sig_x = ta.array([[0, 1], [1, 0]])
sig_y = ta.array([[0, -1j], [1j, 0]])
sig_z = ta.array([[1, 0], [0, -1]])

def make_system(a=1, W=7, L=20, barrier=.5, barrierpos=(0, 3, 17, 20),
                mu=0.3, Delta=0.1, t=1.0):
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.cubic(norbs=2)
    
    def shapeWire(pos):
        x, y, z = pos
        return 0 <= y < 20 and 0 <= x**2 + z**2 < 7
    def shapeLead(pos):
        x, y, z = pos
        return 0 <= y < 1 and 0 <= x**2 + z**2 < 7

    syst = kwant.Builder()
    syst[lat.shape(shapeWire, (0, 0, 0))] = \
        (4 * t - mu) * sig_z + Delta * sig_x

    # Two tunnel barriers
#    syst[(lat(x, y) for x in range(barrierpos[0], barrierpos[1])
#         for y in range(W))] = (4 * t + barrier - mu) * sig_z    
#    syst[(lat(x, y) for x in range(barrierpos[2], barrierpos[3])
#         for y in range(W))] = (4 * t + barrier - mu) * sig_z

    # Hoppings
    syst[lat.neighbors()] = -t * sig_z
    
    # Normal metal lead
    sym_left = kwant.TranslationalSymmetry((0, -a, 0))
    lead = kwant.Builder(sym_left, conservation_law=-sig_z, particle_hole=sig_y)
    lead[lat.shape(shapeLead, (0, 0, 0))] = (4 * t - mu) * sig_z
    lead[lat.neighbors()] = -t * sig_z

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def plotConductanceLead(syst, energies):
    # Compute conductance at the first lead
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        # Conductance is N - R_ee + R_he
        data.append(smatrix.submatrix((0, 0), (0, 0)).shape[0] # N
            - smatrix.transmission((0, 0), (0, 0))      # R_ee l1 to l1
            - smatrix.transmission((1, 0), (0, 0))      # R_ee l2 to l1
            + smatrix.transmission((0, 1), (0, 0)))     # R_he
    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()

def plotConductanceSystem(syst, energies):
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        # Conductance is N - R_ee + R_he
        data.append(smatrix.submatrix((0, 0), (0, 0)).shape[0]
            + smatrix.transmission((0, 0), (1, 0))
            + smatrix.transmission((1, 0), (1, 0))
            + smatrix.transmission((1, 0), (1, 1)))# e to h
    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()

def main():
    syst = make_system()
    plt.rcParams["figure.figsize"] = (10,7)
    kwant.plot(syst)
#    kwant.plot(syst, site_size=0.18, site_lw=0.01, hop_lw=0.1)

    # Finalize the system.
#    syst = syst.finalized()

    # Compute and plot the conductance
#    plotConductanceLead(syst, energies=[0.005 * i for i in range(-80, 80)])
#    plotConductanceSystem(syst, energies=[0.05 * i for i in range(-30, 30)])

if __name__ == '__main__':
    main()