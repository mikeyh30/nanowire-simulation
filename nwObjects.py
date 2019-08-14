#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 01:42:41 2019

@author: Domi
"""
from tqdm import tqdm
import kwant
import tinyarray as ta
import numpy as np
import scipy.sparse.linalg

s0 = np.identity(2)
sZ = np.array([[1., 0.], [0., -1.]])
sX = np.array([[0., 1.], [1., 0.]])
sY = np.array([[0., -1j], [1j, 0.]])

tauZ = ta.array(np.kron(sZ, s0))
tauX = ta.array(np.kron(sX, s0))
tauY = ta.array(np.kron(sY, s0))
sigZ = ta.array(np.kron(s0, sZ))
sigX = ta.array(np.kron(s0, sX))
sigY = ta.array(np.kron(s0, sY))
tauZsigX = ta.array(np.kron(sZ, sX))
tauZsigY = ta.array(np.kron(sZ, sY))
tauYsigY = ta.array(np.kron(sY, sY))

## Functions ##
def makeNISIN(width=5, length=20, barrierLen=1, 
              periodB=.5, M=0.1,
              isWhole=True
              ):
    ## Define site Hopping and functions ##
    def hopX(site0, site1, t, alpha):
        return -t * tauZ + .5j * alpha * tauZsigY
    def hopY(site0, site1, t, alpha):
        return -t * tauZ - .5j * alpha * tauZsigX
        
    def sinuB(theta, M=0.05):
        return M*(sigX*np.cos(theta) - sigY*np.sin(theta))
   
    def onsiteSc(site, muSc, t, B, Delta):
        if periodB == 0:
            return (4 * t - muSc) * tauZ + B * sigX + Delta * tauX
        else:
            theta = 2*np.pi*periodB* \
                (site.pos[0] - barrierLen)/(length - 1 - 2*barrierLen)
            return (4 * t - muSc) * tauZ + B * sigX + Delta * tauX \
                    + sinuB(theta)
    def onsiteNormal(site, mu, t):
        return (4 * t - mu) * tauZ
    def onsiteBarrier(site, mu, t, barrier):
        return (4 * t - mu + barrier) * tauZ
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=4) 
    syst = kwant.Builder()
        
    # S
    syst[(lat(i, j) for i in range(1,length-1) for j in range(width))] \
        = onsiteSc
        
    if isWhole:
        # I's
        syst[(lat(i, j) for i in range(barrierLen) for j in range(width))] = onsiteBarrier
        syst[(lat(i, j) for i in range(length-barrierLen, length)for j in range(width))] = onsiteBarrier
        
        # Hopping:
        syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
        syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    
        # N's
        lead = kwant.Builder(kwant.TranslationalSymmetry((-1,0)),
                             conservation_law=-tauZ,
                             particle_hole=tauYsigY
                             )
        lead[(lat(0, j) for j in range(width))] = onsiteNormal
        lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
        lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
        
        syst.attach_lead(lead)
        syst.attach_lead(lead.reversed())
    else:
        # Hopping:
        syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
        syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    
    return syst.finalized()

## Objects ##
class Nanowire:
    def __init__(self, width=5, length=20, mu=0.0, barrierLen=1, periodB=0, dim=2, M=0.1):
        self.width = width
        self.length = length
        self.mu = mu
        self.barrierLen = barrierLen
        self.periodB = periodB
        self.dim = dim
        self.M = M
        
    def spectrum(self, 
                 bValues=np.linspace(0, 1.0, 101)
                 ):        
        syst = makeNISIN(width=self.width, length=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=False,
                         M=self.M)
        energies = []
        critB = 0
        params = dict(muSc=.0, mu=.3, Delta=.1, alpha=.8, t=1., barrier=2.)
        for i in tqdm(range(np.size(bValues)), 
                      desc= "Spec for periodB = %2.1f" %(self.periodB)
                      ):
            b = bValues[i]
            params["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True,  params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)
            if critB==0 and np.abs(eigs[10] - eigs[9])/2 < 1e-3:
                critB = b
            
        outcome = dict(B=bValues, E=energies, CritB=critB)
        return outcome
    
    def conductances(self, 
                     bValues=np.linspace(0, 1.0, 101),
                     energies=[0.001 * i for i in range(-130, 130)]
                     ):
        syst = makeNISIN(width=self.width, length=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=True,
                         M=self.M)
        data = []
        critB = 0
        params = dict(muSc=.0, mu=.3, Delta=.1, alpha=.8, t=1., barrier=2.)
        for i in tqdm(range(np.size(energies)), 
                      desc= "Cond for periodB = %2.1f" %(self.periodB)
                      ):
            cond = []
            energy = energies[i]
            for b in bValues:
                params["B"] = b
                smatrix = kwant.smatrix(syst, energy, params=params)
                conduct = (
                        smatrix.submatrix((0, 0), (0, 0)).shape[0]  # N
                        - smatrix.transmission((0, 0), (0, 0))      # R_ee
                        + smatrix.transmission((0, 1), (0, 0)))     # R_he
                cond.append(conduct)     # R_he
                if energy == 0 and critB == 0 and np.abs(2 - conduct) < 0.01:
                    critB = b
            data.append(cond)
            
        outcome = dict(B=bValues, BiasV=energies, Cond=data, CritB=critB)
        return outcome
    
    def phaseTransition(self,
                        bValues=np.linspace(0, .5, 101),
                        muValues=np.linspace(0, .5, 101)
                     ):
        syst = makeNISIN(width=self.width, length=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=True,
                         M=self.M)
        criticalPoints = muValues
        params = dict(mu=.3, Delta=.1, alpha=.8, t=1.0, barrier=2.)
        for i in tqdm(range(np.size(muValues)),
                      desc= "Spec for periodB = %2.1f" %(self.periodB)
                          ):
            critB = 0
            params["muSc"] = muValues[i]
            for b in bValues:
                params["B"] = b
                H = syst.hamiltonian_submatrix(sparse=True,  params=params)
                H = H.tocsc()
                eigs = scipy.sparse.linalg.eigsh(H, k=2, sigma=0)
                eigs = np.sort(eigs[0])
                if critB==0 and np.abs(eigs[0] - eigs[1])/2 < 1e-3:
                    critB = b
                    break
            else:
                continue
            criticalPoints[i] = critB
            
        outcome = dict(MuSc=muValues, CritB=criticalPoints)
        return outcome

def main():
    pass

if __name__ == '__main__':
    main()