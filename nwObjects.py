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
def makeNISIN(W=5, L=20, barrierLen=1, periodB=.5, isWhole=True, dim=2):
    ## Define site Hopping and functions ##
    def hopX(site0, site1, t, alpha):
        return -t * tauZ + .5j * alpha * tauZsigY
    def sinuB(theta, B=0.1):
        return B*(sigX*np.cos(theta) - sigY*np.sin(theta))
    if dim == 1:
        def onsiteSc(site, t, B, Delta):
            if periodB == 0:
                return (2 * t) * tauZ + B * sigZ + Delta * tauX
            else:
                theta = 2*np.pi*periodB*(site.pos[0] - barrierLen)/(L - 1 - 2*barrierLen)
                return (2 * t) * tauZ + B * sigZ + Delta * tauX + sinuB(theta)
        def onsiteNormal(site, mu, t):
            return (2 * t - mu) * tauZ
        def onsiteBarrier(site, mu, t, barrier):
            return (2 * t - mu + barrier) * tauZ
        
        ## Create the system and leads ##
        # On each site, electron and hole orbitals.
        lat = kwant.lattice.chain(norbs=4) 
        syst = kwant.Builder()
            
        # S
        syst[(lat(i) for i in range(1,L-1))] = onsiteSc
            
        if isWhole:
            # I's
            syst[(lat(i) for i in range(barrierLen))] = onsiteBarrier
            syst[(lat(i) for i in range(L-barrierLen, L))] = onsiteBarrier
            
            # Hopping:
            syst[lat.neighbors()] = hopX
        
            # N's
            lead = kwant.Builder(kwant.TranslationalSymmetry(*lat.prim_vecs),
                                 conservation_law=-tauZ,
                                 particle_hole=tauYsigY
                                 )
            lead[(lat(0))] = onsiteNormal
            lead[lat.neighbors()] = hopX
            
            syst.attach_lead(lead)
            syst.attach_lead(lead.reversed())
        else:
            # Hopping:
            syst[lat.neighbors()] = hopX
    elif (dim == 2 or dim == 3):
        def hopY(site0, site1, t, alpha):
            return -t * tauZ - .5j * alpha * tauZsigX
        if dim == 2:
            def onsiteSc(site, t, B, Delta):
                if periodB == 0:
                    return (4 * t) * tauZ + B * sigZ + Delta * tauX
                else:
                    theta = 2*np.pi*periodB*(site.pos[0] - barrierLen)/(L - 1 - 2*barrierLen)
                    return (4 * t) * tauZ + B * sigZ + Delta * tauX + sinuB(theta)
            def onsiteNormal(site, mu, t):
                return (4 * t - mu) * tauZ
            def onsiteBarrier(site, mu, t, barrier):
                return (4 * t - mu + barrier) * tauZ
            
            ## Create the system and leads ##
            # On each site, electron and hole orbitals.
            lat = kwant.lattice.square(norbs=4) 
            syst = kwant.Builder()
                
            # S
            syst[(lat(i, j) for i in range(1,L-1) for j in range(W))] = onsiteSc
                
            if isWhole:
                # I's
                syst[(lat(i, j) for i in range(barrierLen) for j in range(W))] = onsiteBarrier
                syst[(lat(i, j) for i in range(L-barrierLen, L)for j in range(W))] = onsiteBarrier
                
                # Hopping:
                syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
                syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
            
                # N's
                lead = kwant.Builder(kwant.TranslationalSymmetry((-1,0)),
                                     conservation_law=-tauZ,
                                     particle_hole=tauYsigY
                                     )
                lead[(lat(0, j) for j in range(W))] = onsiteNormal
                lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
                lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
                
                syst.attach_lead(lead)
                syst.attach_lead(lead.reversed())
            else:
                # Hopping:
                syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
                syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
                
        else:
            def hopZ(site0, site1, t, alpha):
                return -t * tauZ
            def onsiteSc(site, t, B, Delta):
                if periodB == 0:
                    return (6 * t) * tauZ + B * sigX + Delta * tauX
                else:
                    theta = 2*np.pi*periodB*(site.pos[0] - barrierLen)/(L - 1 - 2*barrierLen)
                    return (6 * t) * tauZ + B * sigZ + Delta * tauX + sinuB(theta)
            def onsiteNormal(site, mu, t):
                return (6 * t - mu) * tauZ
            def onsiteBarrier(site, mu, t, barrier):
                return (6 * t - mu + barrier) * tauZ
            
            def wireShape(pos):
                x, y, z = pos
                return 1 <= x < L-1 and 0 <= y**2 + z**2 < W
            def barrier0Shape(pos):
                x, y, z = pos
                return x == 0 and 0 <= y**2 + z**2 < W
            def barrier1Shape(pos):
                x, y, z = pos
                return x == L-1 and 0 <= y**2 + z**2 < W
            
             ## Create the system and leads ##
            # On each site, electron and hole orbitals.
            lat = kwant.lattice.cubic(norbs=4) 
            syst = kwant.Builder()
            
            # S
            syst[lat.shape(wireShape, (0, 0, 0))] = onsiteSc
            
            if isWhole:
                # I's
                syst[lat.shape(barrier0Shape, (0, 0, 0))] = onsiteBarrier
                syst[lat.shape(barrier1Shape, (L-1, 0, 0))] = onsiteBarrier
                
                # Hopping:
                syst[kwant.builder.HoppingKind((1, 0, 0), lat, lat)] = hopX
                syst[kwant.builder.HoppingKind((0, 1, 0), lat, lat)] = hopY
                syst[kwant.builder.HoppingKind((0, 0, 1), lat, lat)] = hopZ
            
                # N's
                lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0, 0)),
                                     conservation_law=-tauZ,
                                     particle_hole=tauYsigY
                                     )
                lead[lat.shape(barrier0Shape, (0, 0, 0))] = onsiteNormal
                lead[kwant.builder.HoppingKind((1, 0, 0), lat, lat)] = hopX
                lead[kwant.builder.HoppingKind((0, 1, 0), lat, lat)] = hopY
                lead[kwant.builder.HoppingKind((0, 0, 1), lat, lat)] = hopZ
            
                syst.attach_lead(lead)
                syst.attach_lead(lead.reversed())
            else:
                syst[kwant.builder.HoppingKind((1, 0, 0), lat, lat)] = hopX
                syst[kwant.builder.HoppingKind((0, 1, 0), lat, lat)] = hopY
                syst[kwant.builder.HoppingKind((0, 0, 1), lat, lat)] = hopZ
    else:
         print("Wrong dimensions")
         pass
    
    return syst.finalized()

## Objects ##
class Nanowire:
    def __init__(self, width=5, length=20, barrierLen=1, periodB=0, dim=2):
        self.width = width
        self.length = length
        self.barrierLen = barrierLen
        self.periodB = periodB
        self.dim = dim
        
    def spectrum(self, 
                 bValues=np.linspace(0, 1.0, 101)
                 ):        
        syst = makeNISIN(W=self.width, L=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=False,
                         dim=self.dim)
        energies = []
        params = dict(
                mu=.3, Delta=.1, alpha=.8, t=1.0, barrier=2.
                )
        for i in tqdm(range(np.size(bValues)), 
                      desc= "Spec for periodB = %2.1f" %(self.periodB)
                      ):
            params["B"] = bValues[i]
            H = syst.hamiltonian_submatrix(sparse=True,  params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            energies.append(np.sort(eigs[0]))
            
        outcome = dict(B=bValues, E=energies)
        return outcome
    
    def conductances(self, 
                     bValues=np.linspace(0, 1.0, 101),
                     energies=[0.001 * i for i in range(-130, 130)]
                     ):
        syst = makeNISIN(W=self.width, L=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=True,
                         dim=self.dim)
        data = []
        critB = 0
        params = dict(
                mu=.3, Delta=.1, alpha=.8, t=1.0, barrier=2.
                )
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
                if energy == 0 and critB == 0:
                    if np.abs(2 - conduct) < 0.01:
                        critB = b
            data.append(cond)
            
        outcome = dict(B=bValues, BiasV=energies, Cond=data, CritB=critB)
        return outcome
    
    def criticalValue(self, 
                     bValues=np.linspace(0, 1., 101)
                     ):
        syst = makeNISIN(W=self.width, L=self.length, 
                         barrierLen=self.barrierLen, 
                         periodB=self.periodB, 
                         isWhole=True,
                         dim=self.dim)
        params = dict(
                mu=.3, Delta=.1, alpha=.8, t=1.0, barrier=2.
                )
        
        for i in tqdm(range(np.size(bValues)), 
                      desc= "Spec for periodB = %2.1f" %(self.periodB)
                      ):
            critB = bValues[i]
            params["B"] = critB
            smatrix = kwant.smatrix(syst, 0, params=params)
            cond = (smatrix.submatrix((0, 0), (0, 0)).shape[0]  # N
                    - smatrix.transmission((0, 0), (0, 0))      # R_ee
                    + smatrix.transmission((0, 1), (0, 0)))     # R_he
            if np.abs(2 - cond) < 0.01:
                return critB
            
        return 0

def main():
    pass

if __name__ == '__main__':
    main()