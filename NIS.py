#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 01:33:22 2019

@author: Domi
"""

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

def makeNIS1D(
        mu=.5, B=0., Delta=.1, alpha=0.15, t=1.0,
        barrier = 1.
        ):
    onsiteSc      = (2 * t) * tau_z + B * sig_z + Delta * tau_x
    onsiteNormal  = (2 * t - mu) * tau_z + B * sig_z
    onsiteBarrier = (2 * t - mu + barrier) * tau_z
    hop           = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.chain(norbs=4) 
    syst = kwant.Builder()
        
    # Barriers
    syst[lat(0)] = onsiteBarrier
    
    # Normal metal lead
    lead0 = kwant.Builder(kwant.TranslationalSymmetry(*lat.prim_vecs), 
                          conservation_law=-tau_z
                          )
    lead0[lat(0)] = onsiteNormal
    lead0[lat.neighbors()] = hop
    
    lead1 = kwant.Builder(kwant.TranslationalSymmetry(*lat.prim_vecs))
    lead1[lat(0)] = onsiteSc
    lead1[lat.neighbors()] = hop

    syst.attach_lead(lead0.reversed())
    syst.attach_lead(lead1)

    return syst

def makeNIS2D(
        W=5,
        mu=.5, B=0., Delta=.1, alpha=.8, t=1.0,
        barrier = 1.
        ):
    onsiteSc      = (4 * t) * tau_z + B * sig_z + Delta * tau_x
    onsiteNormal  = (4 * t - mu) * tau_z + B * sig_z
    onsiteBarrier = (4 * t - mu + barrier) * tau_z
    hop           = -t * tau_z - .5j * alpha * tau_zsig_x
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=4) 
    syst = kwant.Builder()
        
    # Barriers
    syst[(lat(0, j) for j in range(W))] = onsiteBarrier
    syst[lat.neighbors()] = hop
    
    # Normal metal lead
    lead0 = kwant.Builder(kwant.TranslationalSymmetry((-1,0)), 
                          conservation_law=-tau_z
                          )
    lead0[(lat(0, j) for j in range(W))] = onsiteNormal
    lead0[lat.neighbors()] = hop
    
    lead1 = kwant.Builder(kwant.TranslationalSymmetry((1,0)))
    lead1[(lat(0, j) for j in range(W))] = onsiteSc
    lead1[lat.neighbors()] = hop

    syst.attach_lead(lead0)
    syst.attach_lead(lead1)

    return syst

def plotConductance(syst, energies):
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

def main():
#    syst = makeNIS1D()
    syst = makeNIS2D()
    plt.rcParams["figure.figsize"] = (8,5)
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # Compute and plot the conductance
    plotConductance(syst, energies=[0.005 * i for i in range(-50, 200)])

if __name__ == '__main__':
    main()