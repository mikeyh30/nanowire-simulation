#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 03:32:47 2019

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

def makeNISIN1D(
        L=20, # B = 0.22
        mu=.15, B=0., Delta=.1, alpha=.1, t=1.0,
        barrier = .5
        ):
    onsiteSc      = (2 * t) * tau_z + B * sig_z + Delta * tau_x
    onsiteNormal  = (2 * t - mu) * tau_z
    onsiteBarrier = (2 * t - mu + barrier) * tau_z
    hop           = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.chain(norbs=4) 
    syst = kwant.Builder()
    
    # S
    syst[(lat(i) for i in range(1,L-1))] = onsiteSc
    
    # I's
    syst[lat(0)] = onsiteBarrier
    syst[lat(L-1)] = onsiteBarrier
    syst[lat.neighbors()] = hop
    
    # N's
    lead = kwant.Builder(kwant.TranslationalSymmetry(*lat.prim_vecs),
                         conservation_law=-tau_z
                         )
    lead[lat(0)] = onsiteNormal
    lead[lat.neighbors()] = hop

    syst.attach_lead(lead.reversed())
    syst.attach_lead(lead)

    return syst

def makeNISIN2D(
        W=5,
        L=20,
        mu=.3, B=0., Delta=.1, alpha=.8, t=1.0,
        barrier = .5
        ):
    onsiteSc      = (4 * t) * tau_z + B * sig_z + Delta * tau_x
    onsiteNormal  = (4 * t - mu) * tau_z
    onsiteBarrier = (4 * t - mu + barrier) * tau_z
    hop           = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=4) 
    syst = kwant.Builder()
    
    # S
    syst[(lat(i, j) for i in range(1,L-1) for j in range(W))] = onsiteSc
    
    barrierLen = 1;
    # I's
    syst[(lat(i, j) for i in range(barrierLen) for j in range(W))] = onsiteBarrier
    syst[(lat(i, j) for i in range(L-barrierLen, L)for j in range(W))] = onsiteBarrier
    syst[lat.neighbors()] = hop
    
    # N's
    lead = kwant.Builder(kwant.TranslationalSymmetry((-1,0)),
                         conservation_law=-tau_z
                         )
    lead[(lat(0, j) for j in range(W))] = onsiteNormal
    lead[lat.neighbors()] = hop

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def plotConductanceSc(syst, energies):
    # Compute conductance at the first lead
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        # Conductance is N - R_ee + R_he
        data.append(
                smatrix.submatrix((0, 0), (0, 0)).shape[0]
                    )     # R_he
    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()
    
def plotConductance(syst, energies):
    # Compute conductance at the first lead
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        # Conductance is N - R_ee + R_he
        data.append(
#                smatrix.submatrix((0, 0), (0, 0)).shape[0]
#                - smatrix.transmission((0, 0), (0, 0))
                + smatrix.transmission((1, 1), (0, 0))
                    )     # R_he
    plt.figure()
    plt.plot(energies, data)
    plt.xlabel("energy [t]")
    plt.ylabel("conductance [e^2/h]")
    plt.show()

def main():
#    syst = makeNISIN1D(B=0.)
    syst = makeNISIN2D(L=20, B=0.2) #.25
    plt.rcParams["figure.figsize"] = (8,5)
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # Compute and plot the conductance
#    plotConductanceSc(syst, energies=[0.01 * i for i in range(-50, 50)])
    plotConductance(syst, energies=[0.001 * i for i in range(-160, 160)])

if __name__ == '__main__':
    main()