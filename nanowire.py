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

def make_system(W=7, L=20, barrier=5, barrierPos=(0, 3, 17, 20),
                mu=1.9, Delta=.2, t=1.0):
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.cubic(norbs=2)
    
    def shapeWire(pos):
        x, y, z = pos
        return 0 <= y < L and x**2 + z**2 < W
    def shapeBarrier0(pos):
        x, y, z = pos
        return barrierPos[0] <= y < barrierPos[1] and x**2 + z**2 < 7
    def shapeBarrier1(pos):
        x, y, z = pos
        return barrierPos[2] <= y < barrierPos[3] and x**2 + z**2 < 7
    def shapeLead(pos):
        x, y, z = pos
        return 0 <= y < 1 and x**2 + z**2 < 7
    
    # Initialise empty tight-binding model
    syst = kwant.Builder()
    
    # Super
    syst[lat.shape(shapeWire, (0, 0, 0))] = \
        (6 * t - mu) * sig_z + Delta * sig_x

    # Two tunnel barriers
    syst[lat.shape(shapeBarrier0, (0, 0, 0))] = (6 * t + barrier - mu) * sig_z
    syst[lat.shape(shapeBarrier1, (0, L-1, 0))] = (6 * t + barrier - mu) * sig_z

    # Hoppings
    syst[lat.neighbors()] = -t * sig_z
    
    # Normal metal lead
    sym_left = kwant.TranslationalSymmetry((0, -1, 0))
    lead = kwant.Builder(sym_left, conservation_law=-sig_z, particle_hole=sig_y)
    lead[lat.shape(shapeLead, (0, 0, 0))] = (6 * t - mu) * sig_z
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

def family_color(site):
    x, y, z = site.pos
    if 0 <= y < 3 or 17 <=y < 20:
        return 'green'
    else:
        return 'blue'

def main():
    syst = make_system()
    plt.rcParams["figure.figsize"] = (12,8)
    kwant.plot(syst, site_color=family_color)
#    kwant.plot(syst, site_size=0.18, site_lw=0.01, hop_lw=0.1)

    # Finalize the system.
    syst = syst.finalized()

    # Compute and plot the conductance
    plotConductanceLead(syst, energies=[0.005 * i for i in range(-120, 120)])
    plotConductanceSystem(syst, energies=[0.05 * i for i in range(-50, 50)])

if __name__ == '__main__':
    main()