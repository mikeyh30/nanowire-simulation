#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 13:59:55 2019

@author: Domi
"""
import sys
sys.path.append("/object")
from nanowireObject import *

import kwant
import numpy as np
import scipy.sparse.linalg
import matplotlib.pyplot as plt

import tinyarray as ta

class SimpleNamespace(object):
    """A simple container for parameters."""
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

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

def onsite(site, p):
    return tau_z * (p.mu - 2 * p.t) + \
        sig_z * p.B + tau_x * p.Delta
        
def hopping(site0, site1, p):
    return tau_z * p.t + 1j * tau_zsig_x * p.alpha

def make_system(l=70):
    sys = kwant.Builder()
    lat = kwant.lattice.chain()
    sys[(lat(x) for x in range(l))] = onsite
    sys[lat.neighbors()] = hopping
    return sys.finalized()

sys = make_system()

# Calculate and plot lowest eigenenergies in B-field.
B_values = np.linspace(0, 0.6, 80)
energies = []
params = SimpleNamespace(
    t=1, mu=-0.1, alpha=0.05, Delta=0.2)
for params.B in B_values:
    H = sys.hamiltonian_submatrix(
        args=[params], sparse=True)
    H = H.tocsc()
    eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
    energies.append(np.sort(eigs[0]))
plt.plot(B_values, energies)
plt.show()