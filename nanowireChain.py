#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:59:55 2019

@author: Domi
"""

import numpy as np
import tinyarray as ta
import scipy.sparse.linalg

import kwant
import matplotlib.pyplot as plt

s_0 = np.identity(2)
s_z = np.array([[1., 0.], [0., -1.]])
s_x = np.array([[0., 1.], [1., 0.]])
s_y = np.array([[0., -1j], [1j, 0.]])

tau_z = ta.array(np.kron(s_z, s_0))
tau_x = ta.array(np.kron(s_x, s_0))
tau_y = ta.array(np.kron(s_y, s_0))
sig_z = ta.array(np.kron(s_0, s_z))
sig_x = ta.array(np.kron(s_0, s_x))
sig_y = ta.array(np.kron(s_0, s_y))
tau_zsig_x = ta.array(np.kron(s_z, s_x))

def onsite(site, mu, t, B, Delta):
    return tau_z * (mu - 2 * t) + \
        sig_z * B + tau_x * Delta
        
def hopping(site0, site1, t, alpha):
    return tau_z * t + 1j * tau_zsig_x * alpha

def make_system(l=70):
    sys = kwant.Builder()
    lat = kwant.lattice.chain()
    sys[(lat(x) for x in range(l))] = onsite
    sys[lat.neighbors()] = hopping
    kwant.plot(sys)
    return sys.finalized()

sys = make_system()

# Calculate and plot lowest eigenenergies in B-field.
B_values = np.linspace(0, 0.6, 80)
energies = []
params = dict(
        t=1, mu=-0.1, alpha=0.05, Delta=0.2
        )
for B in B_values:
    params["B"] = B
    H = sys.hamiltonian_submatrix(sparse=True,  params=params)
#    print(tupleH[2][1]) # to check use: return_norb=True,
    H = H.tocsc()
    eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0) # near zero!! this is zero modes
    energies.append(np.sort(eigs[0]))
plt.plot(B_values, energies)
plt.xlabel("B")
plt.ylabel("energies")
plt.show()