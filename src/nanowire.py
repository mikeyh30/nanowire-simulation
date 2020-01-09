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
from nanomagnet_field import rick_fourier
from scipy.constants import physical_constants, hbar
from transport_model import NISIN

bohr_magneton = 1  # physical_constants['Bohr magneton'][0]
lattice_constant_InAs = 1  # 6.0583E-10 # might need to change this.

## Objects ##
class Nanowire:
    def __init__(
        self,
        width=5,
        noMagnets=5,
        dim=2,
        barrier_length=1,
        effective_mass=0.5,
        M=0.05,
        muSc=0.0,
        alpha=0.8,
        addedSinu=False,
        stagger_ratio=0.5,
        mu=0.3,
        delta=0.1,
        barrier=2.0,
    ):
        # Wire Physical Properties
        self.width = width
        self.noMagnets = noMagnets
        self.dim = dim
        self.barrier_length = barrier_length

        # Superconducting components
        self.t = 0.5 / effective_mass # (hbar**2)/(2*effective_mass*lattice_constant_InAs)
        self.M = M
        self.muSc = muSc
        self.alpha = alpha
        self.addedSinu = addedSinu

        # Nanomagnet properties
        self.stagger_ratio = stagger_ratio

        # Previously hard-coded parameters
        self.mu = mu  # how is this different from muSc?
        self.delta = delta
        self.barrier = barrier

        # System
        """self.system = NISIN(width=self.width, noMagnets=self.noMagnets, 
                                barrier_length=self.barrier_length, M=self.M,
                                addedSinu=self.addedSinu, 
                                stagger_ratio=self.stagger_ratio
                                )
        """

    def spectrum(self, bValues=np.linspace(0, 1.0, 201)):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
        )
        energies = []
        critB = 0
        params = dict(
            muSc=self.muSc,
            mu=self.mu,
            Delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            addedSinu=self.addedSinu,
            M=self.M,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length
        )
        for i in tqdm(
            range(np.size(bValues)),
            desc="Spec: NoMagnets = %i, added? %r" % (self.noMagnets, self.addedSinu),
        ):
            b = bValues[i]
            params["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True, params=params)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)
            if critB == 0 and np.abs(eigs[10] - eigs[9]) / 2 < 1e-3:
                critB = b

        outcome = dict(B=bValues, E=energies, CritB=critB)
        return outcome

    def conductances(
        self,
        bValues=np.linspace(0, 1.0, 201),
        energies=[0.001 * i for i in range(-120, 120)],
    ):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length
        )
        data = []
        critB = 0
        params = dict(
            muSc=self.muSc,
            mu=self.mu,
            Delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            M=self.M,
            addedSinu=self.addedSinu,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length
        )
        for i in tqdm(
            range(np.size(energies)),
            desc="Cond: NoMagnets = %i, added? %r" % (self.noMagnets, self.addedSinu),
        ):
            cond = []
            energy = energies[i]
            for b in bValues:
                params["B"] = b
                smatrix = kwant.smatrix(syst, energy, params=params)
                conduct = (
                    smatrix.submatrix((0, 0), (0, 0)).shape[0]  # N
                    - smatrix.transmission((0, 0), (0, 0))  # R_ee
                    + smatrix.transmission((0, 1), (0, 0))
                )  # R_he
                cond.append(conduct)
                if (
                    np.isclose(energy, 0, rtol=1e-6)
                    and critB == 0
                    and np.abs(2 - conduct) < 0.01
                ):
                    critB = b
            data.append(cond)

        outcome = dict(B=bValues, BiasV=energies, Cond=data, CritB=critB)
        return outcome

    def criticals(
        self, bValues=np.linspace(0, 0.4, 81), muValues=np.linspace(0, 1.0, 101)
    ):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
        )
        criticalPoints = []
        params = dict(
            mu=self.mu,
            Delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            M=self.M,
            addedSinu=self.addedSinu,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length
        )
        for i in tqdm(
            range(np.size(muValues)),
            desc="Crit: NoMagnets = %i, added? %r" % (self.noMagnets, self.addedSinu),
        ):
            params["muSc"] = muValues[i]
            for b in bValues:
                params["B"] = b
                H = syst.hamiltonian_submatrix(sparse=True, params=params)
                H = H.tocsc()
                eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
                eigs = np.sort(eigs[0])
                if np.abs(eigs[9] - eigs[10]) / 2 < 1e-3:
                    criticalPoints.append(b)
                    break
            else:
                continue

        outcome = dict(MuSc=muValues, CritB=criticalPoints)
        return outcome

    def phaseTransitions(
        self, bValues=np.linspace(0, 0.4, 81), muValues=np.linspace(0, 1.0, 101)
    ):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
        )
        phases = []
        params = dict(
            mu=self.mu,
            Delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            M=self.M,
            addedSinu=self.addedSinu,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length
        )
        for i in tqdm(
            range(np.size(bValues)),
            desc="Phas: NoMagnets = %i, added? %r" % (self.noMagnets, self.addedSinu),
        ):
            params["B"] = bValues[i]
            transitions = []
            for mu in muValues:
                params["muSc"] = mu
                H = syst.hamiltonian_submatrix(sparse=True, params=params)
                H = H.tocsc()
                eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
                eigs = np.sort(eigs[0])
                if 0.5 * np.abs(eigs[9] - eigs[10]) < 1e-3:
                    transitions.append(1)
                else:
                    transitions.append(0)

            phases.append(transitions)

        outcome = dict(B=bValues, MuSc=muValues, Phases=phases)
        return outcome

    def phaseAid(
        self, bValues=np.linspace(0.0, 0.3, 61), muValues=np.linspace(0.2, 0.5, 61)
    ):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
        )
        energies0 = []
        energies1 = []
        params = dict(
            muSc=muValues[0],
            mu=self.mu,
            Delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            M=self.M,
            addedSinu=self.addedSinu,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length
        )
        for i in tqdm(
            range(np.size(bValues)), desc="Number of magnets = %i" % (self.noMagnets)
        ):
            b = bValues[i]
            params["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True, params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            energies0.append(np.sort(eigs[0]))

        for i in tqdm(
            range(np.size(muValues)), desc="Number of magnets = %i" % (self.noMagnets)
        ):
            mu = muValues[i]
            params["muSc"] = mu
            H = syst.hamiltonian_submatrix(sparse=True, params=params)
            H = H.tocsc()
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            energies1.append(np.sort(eigs[0]))

        outcome = dict(B=bValues, MuSc=muValues, Eb=energies0, Em=energies1)
        return outcome

    def plot(self):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
        )

        return kwant.plotter.plot(syst, show=False)


def main():
    pass


if __name__ == "__main__":
    main()
