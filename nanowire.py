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
from nanomagnet_field import rick_sinusiod
from scipy.constants import physical_constants

bohr_magneton = physical_constants['Bohr magneton'][0]

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
def makeNISIN(width=7, noMagnets=5, barrierLen=1, M=0.05,
              addedSinu=False, isWhole=True, stagger_ratio=0.5
              ):
    length = 8*noMagnets - 2 + 2*barrierLen # "2*(noMagnets - 1)"

    ## Define site Hopping and functions ##
    def hopX(site0, site1, t, alpha):
        return -t * tauZ + 1j * alpha * tauZsigY
    def hopY(site0, site1, t, alpha):
        return -t * tauZ - 1j * alpha * tauZsigX
        
    def sinuB(theta,stagger_ratio):
        ssin, scos = rick_sinusiod(theta)
        return sigY*scos + sigX*ssin
    
    # This is the onsite Hamiltonian, this is where the B-field can be varied.
    def onsiteSc(site, muSc, t, B, Delta, gfactor=10):
        if addedSinu:
            counter = np.mod(site.pos[0]-1-barrierLen, 16)
            if -1 < counter < 4:
                theta = 0
            elif 3 < counter < 8:
                theta = .2*(counter - 3)*np.pi
            elif 7 < counter < 12:
                theta = np.pi
            else:
                theta = .2*(counter - 6)*np.pi
                    
            return (4 * t - muSc) * tauZ + 0.5*gfactor*bohr_magneton*B*sigX + Delta * tauX \
                    + 0.5*gfactor*bohr_magneton*M*sinuB(theta,stagger_ratio)
        else:
            return (4 * t - muSc) * tauZ + 0.5*gfactor*bohr_magneton*B*sigX + Delta * tauX
    def onsiteNormal(site, mu, t):
        return (4 * t - mu) * tauZ
    def onsiteBarrier(site, mu, t, barrier):
        return (4 * t - mu + barrier) * tauZ
    
    # On each site, electron and hole orbitals.
    lat = kwant.lattice.square(norbs=4) 
    syst = kwant.Builder()
        
    # S
    syst[(lat(i, j) 
            for i in range(barrierLen,length-barrierLen) 
            for j in range(width)) ] = onsiteSc
        
    if isWhole:
        # I's
        syst[(lat(i, j) 
                for i in range(barrierLen) 
                for j in range(width)) ] = onsiteBarrier
        syst[(lat(i, j) 
                for i in range(length-barrierLen, length) 
                for j in range(width)) ] = onsiteBarrier
        
        # Hopping:
        syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
        syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    
        # N's
        lead = kwant.Builder(kwant.TranslationalSymmetry((-1,0)),
                             conservation_law=-tauZ,
                             particle_hole=tauYsigY )
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
    def __init__(self, width=5, noMagnets=5, dim=2, barrierLen=1, 
                 effective_mass=.5, M=0.05, muSc=.0, alpha=.8, addedSinu=False,
                 stagger_ratio=0.5, mu=0.3, delta=0.1, barrier=2.0
                 ):
        # Wire Physical Properties
        self.width=width
        self.noMagnets=noMagnets
        self.dim = dim
        self.barrierLen=barrierLen
        
        # Superconducting components
        self.t=.5/effective_mass # check that this definition is correct.
        self.M=M
        self.muSc=muSc
        self.alpha=alpha
        self.addedSinu = addedSinu

        # Nanomagnet properties
        self.stagger_ratio = stagger_ratio

        # Previously hard-coded parameters
        self.mu = mu # how is this different from muSc?
        self.delta = delta
        self.barrier = barrier

        # System
        '''self.system = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                                barrierLen=self.barrierLen, M=self.M,
                                addedSinu=self.addedSinu, isWhole=False,
                                stagger_ratio=self.stagger_ratio
                                )
        '''

    def spectrum(self, bValues=np.linspace(0, 1.0, 201)):        
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=False,
                         stagger_ratio=self.stagger_ratio
                         )
        energies = []
        critB = 0
        params = dict(muSc=self.muSc, mu=self.mu, Delta=self.delta, alpha=self.alpha, 
                      t=self.t, barrier=self.barrier)
        for i in tqdm(range(np.size(bValues)), 
                      desc="Spec: NoMagnets = %i, added? %r" 
                      %(self.noMagnets, self.addedSinu)
                      ):
            b = bValues[i]
            params["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True,  params=params)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)
            if critB==0 and np.abs(eigs[10] - eigs[9])/2 < 1e-3:
                critB = b
            
        outcome = dict(B=bValues, E=energies, CritB=critB)
        return outcome
    
    def conductances(self, 
                     bValues=np.linspace(0, 1.0, 201),
                     energies=[0.001 * i for i in range(-120, 120)]
                     ):
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=True,
                         stagger_ratio=self.stagger_ratio
                         )
        data = []
        critB = 0
        params = dict(muSc=self.muSc, mu=self.mu, Delta=self.delta, alpha=self.alpha, 
                      t=self.t, barrier=self.barrier)
        for i in tqdm(range(np.size(energies)), 
                      desc="Cond: NoMagnets = %i, added? %r" 
                      %(self.noMagnets, self.addedSinu)
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
                cond.append(conduct)
                if np.isclose(energy,0,rel_tol=1E-6) and critB == 0 and np.abs(2 - conduct) < 0.01:
                    critB = b
            data.append(cond)
            
        outcome = dict(B=bValues, BiasV=energies, Cond=data, CritB=critB)
        return outcome
    
    def criticals(self,
                        bValues=np.linspace(0, .4, 81),
                        muValues=np.linspace(0, 1.0, 101)
                     ):
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=False,
                         stagger_ratio=self.stagger_ratio
                         )
        criticalPoints = []
        params = dict(mu=self.mu, Delta=self.delta, alpha=self.alpha, 
                      t=self.t, barrier=self.barrier)
        for i in tqdm(range(np.size(muValues)),
                      desc="Crit: NoMagnets = %i, added? %r" 
                      %(self.noMagnets, self.addedSinu)
                      ):
            params["muSc"] = muValues[i]
            for b in bValues:
                params["B"] = b
                H = syst.hamiltonian_submatrix(sparse=True,  params=params)
                H = H.tocsc()
                eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
                eigs = np.sort(eigs[0])
                if np.abs(eigs[9] - eigs[10])/2 < 1e-3:
                    criticalPoints.append(b)
                    break
            else:
                continue
            
        outcome = dict(MuSc=muValues, CritB=criticalPoints)
        return outcome
    
    def phaseTransitions(self,
                        bValues=np.linspace(0, .4, 81),
                        muValues=np.linspace(0, 1.0, 101)
                     ):
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=False,
                         stagger_ratio=self.stagger_ratio
                         )
        phases = []
        params = dict(mu=self.mu, Delta=self.delta, alpha=self.alpha, 
                      t=self.t, barrier=self.barrier)
        for i in tqdm(range(np.size(bValues)),
                      desc="Phas: NoMagnets = %i, added? %r" 
                      %(self.noMagnets, self.addedSinu)
                      ):
            params["B"] = bValues[i]
            transitions = []
            for mu in muValues:
                params["muSc"] = mu
                H = syst.hamiltonian_submatrix(sparse=True,  params=params)
                H = H.tocsc()
                eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
                eigs = np.sort(eigs[0])
                if .5*np.abs(eigs[9] - eigs[10]) < 1e-3:
                    transitions.append(1)
                else:
                    transitions.append(0)
                    
            phases.append(transitions)
            
        outcome = dict(B=bValues, MuSc=muValues, Phases=phases)
        return outcome
    
    def phaseAid(self, 
                 bValues=np.linspace(.0, .3, 61),
                 muValues=np.linspace(.2, .5, 61)
                 ):
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=False,
                         stagger_ratio=self.stagger_ratio
                         )
        energies0 = []
        energies1 = []
        params = dict(muSc=muValues[0], mu=self.mu, Delta=self.delta, alpha=self.alpha, 
                      t=self.t, barrier=self.barrier)
        for i in tqdm(range(np.size(bValues)), 
                      desc="Number of magnets = %i" %(self.noMagnets)
                      ):
            b = bValues[i]
            params["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True,  params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            energies0.append(np.sort(eigs[0]))
            
        for i in tqdm(range(np.size(muValues)), 
                      desc="Number of magnets = %i" %(self.noMagnets)
                      ):
            mu = muValues[i]
            params["muSc"] = mu
            H = syst.hamiltonian_submatrix(sparse=True,  params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            energies1.append(np.sort(eigs[0]))

            
        outcome = dict(B=bValues, MuSc=muValues, Eb=energies0, Em=energies1)
        return outcome

    def plot(self):
        syst = makeNISIN(width=self.width, noMagnets=self.noMagnets, 
                         barrierLen=self.barrierLen, M=self.M,
                         addedSinu=self.addedSinu, isWhole=False,
                         stagger_ratio=self.stagger_ratio
                         )
        
        return kwant.plotter.plot(syst,show=False)

def main():
    pass

if __name__ == '__main__':
    main()