#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 03:16:59 2019

@author: Domi
"""

import kwant
import tinyarray as ta
import numpy as np
import matplotlib.pyplot as plt

s_0 = np.identity(2)
s_z = np.array([[1., 0.], [0., -1.]])
s_x = np.array([[0., 1.], [1., 0.]])
s_y = np.array([[0., -1j], [1j, 0.]])

tau_z = ta.array(np.kron(s_z, s_0))
tau_x = ta.array(np.kron(s_x, s_0))
tau_y = ta.array(np.kron(s_y, s_0))
sig_z = ta.array(np.kron(s_0, s_z))
tau_zsig_x = ta.array(np.kron(s_z, s_x))

def make2DSystem(W=10, L=100, barrierPos=(0, 2, 98, 100),
                mu=1., B=5, Delta=.2, alpha=0.4, t=1.0):
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=4) 
    
    # Initialise empty tight-binding model
    syst = kwant.Builder()
    
    # Super
    syst[(lat(x, y) for x in range(L) for y in range(W))] = \
        (4 * t) * tau_z + B * sig_z + Delta * tau_x
        
    barrier = 1.0 + mu
    # Two tunnel barriers
    syst[(lat(x, y) for x in range(barrierPos[0], barrierPos[1])
         for y in range(W))] = (4 * t + barrier) * tau_z  
    syst[(lat(x, y) for x in range(barrierPos[2], barrierPos[3])
         for y in range(W))] = (4 * t + barrier) * tau_z

    # Hoppings
    syst[lat.neighbors()] = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # Normal metal lead
    sym_left = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_left, conservation_law=-tau_z)
    lead[(lat(0, j) for j in range(W))] = (4 * t - mu) * tau_z
    lead[lat.neighbors()] = -t * tau_z

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def make3DSystem(W=7, L=20, barrierPos=(0, 2, 18, 20),
                mu=1., B=5, Delta=.2, alpha=0.4, t=1.0):
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.cubic(norbs=4)
    
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
        (6 * t - mu) * tau_z + B * sig_z + Delta * tau_x
        
    barrier = .2 + mu
    # Two tunnel barriers
    syst[lat.shape(shapeBarrier0, (0, 0, 0))] = \
        (6 * t + barrier) * tau_z
    syst[lat.shape(shapeBarrier1, (0, L-1, 0))] = \
        (6 * t + barrier) * tau_z

    # Hoppings
    syst[lat.neighbors()] = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # Normal metal lead
    sym_left = kwant.TranslationalSymmetry((0, -1, 0))
    lead = kwant.Builder(sym_left, conservation_law=-tau_z, particle_hole=tau_y)
    lead[lat.shape(shapeLead, (0, 0, 0))] = (6 * t - mu) * tau_z
    lead[lat.neighbors()] = -t * tau_z

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def plotConductanceSystem(syst, energies):
    # Compute conductance at the first lead
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        # Conductance is N - R_ee + R_he
        data.append(smatrix.submatrix((0, 0), (0, 0)).shape[0]
                    - smatrix.transmission((0, 0), (0, 0))
                    + smatrix.transmission((0, 1), (0, 0))
                    )     # R_he
    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()

def family3Dcolor(site):
    x, y, z = site.pos
    if 0 <= y < 3 or 17 <=y < 20:
        return 'green'
    else:
        return 'blue'
    
def family2Dcolor(site):
    x, y = site.pos
    if 0 <= x < 2 or 98 <= x < 100:
        return 'green'
    else:
        return 'blue'

def main():
    syst = make2DSystem(W=8, barrierPos=(0, 2, 98, 100),
        mu=0.5, B=0., Delta=.1, alpha=0.15, t=1.0)
    plt.rcParams["figure.figsize"] = (8,5)
    kwant.plot(syst, site_color=family2Dcolor)
#    kwant.plot(syst, site_size=0.18, site_lw=0.01, hop_lw=0.1)

    # Finalize the system.
    syst = syst.finalized()

    # Compute and plot the conductance
    plotConductanceSystem(syst, energies=[0.005 * i for i in range(-80, 80)])

if __name__ == '__main__':
    main()