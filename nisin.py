#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:35:35 2019

@author: Domi
"""

import kwant
import tinyarray as ta
#import numpy as np
import matplotlib.pyplot as plt

sig_x = ta.array([[0, 1], [1, 0]])
sig_y = ta.array([[0, -1j], [1j, 0]])
sig_z = ta.array([[1, 0], [0, -1]])

def make_system(a=1, W=7, L=21, barrier=.5, barrierpos=(0, 2, 18, 20),
                mu=0.3, Delta=0.1, t=1.0):
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=2)
    # Create an empty tight binding model.
    syst = kwant.Builder()
    
    # Scattering region
    syst[(lat(x, y) for x in range(L) for y in range(W))] = \
        (4 * t - mu) * sig_z + Delta * sig_x

    # Two tunnel barriers
    syst[(lat(x, y) for x in range(barrierpos[0], barrierpos[1])
         for y in range(W))] = (4 * t + barrier - mu) * sig_z    
    syst[(lat(x, y) for x in range(barrierpos[2], barrierpos[3])
         for y in range(W))] = (4 * t + barrier - mu) * sig_z

    # Hoppings
    syst[lat.neighbors()] = -t * sig_z
    
    # Normal metal lead
    sym_left = kwant.TranslationalSymmetry((-a, 0))
    lead = kwant.Builder(sym_left, conservation_law=-sig_z, particle_hole=sig_y)
    lead[(lat(0, j) for j in range(W))] = (4 * t - mu) * sig_z
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
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # Compute and plot the conductance
    plotConductanceLead(syst, energies=[0.005 * i for i in range(-80, 80)])
    plotConductanceSystem(syst, energies=[0.05 * i for i in range(-50, 200)])

if __name__ == '__main__':
    main()