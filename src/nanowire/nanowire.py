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
from nanowire.nanomagnet_field import rick_fourier
from nanowire.transport_model import NISIN, barrier_region


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
        alpha_R=0.8,
        added_sinusoid=False,
        stagger_ratio=0.5,
        mu=0.3,
        delta=0.1,
        barrier=2.0,
        hopping_distance=1,
        bohr_magneton=1,
        gfactor=1,
    ):
        # Wire Physical Properties
        self.width = width
        self.noMagnets = noMagnets
        self.dim = dim
        self.barrier_length = barrier_length

        # Superconducting components
        self.t = 3.83 / (
            effective_mass * (hopping_distance ** 2)
        )  # (hbar**2)/(2*effective_mass*electron_mass*(self.hopping_distance**2))# 0.5 / effective_mass
        self.M = M
        self.muSc = muSc
        self.alpha = alpha_R / hopping_distance
        self.added_sinusoid = added_sinusoid

        # Nanomagnet properties
        self.stagger_ratio = stagger_ratio

        # Previously hard-coded parameters
        self.mu = mu  # how is this different from muSc?
        self.delta = delta
        self.barrier = barrier

        self.gfactor = gfactor
        self.bohr_magneton = bohr_magneton
        self.hopping_distance = hopping_distance

    def spectrum(self, bValues=np.linspace(0, 1.0, 201)):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
            hopping_distance=self.hopping_distance,
        )
        energies = []
        critB = 0
        params = dict(
            muSc=self.muSc,
            mu=self.mu,
            delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            added_sinusoid=self.added_sinusoid,
            M=self.M,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length,
            gfactor=self.gfactor,
            bohr_magneton=self.bohr_magneton,
            hopping_distance=self.hopping_distance,
        )
        for i in tqdm(range(np.size(bValues)), desc="Spec",):
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
            barrier_length=self.barrier_length,
            hopping_distance=self.hopping_distance,
        )
        data = []
        critB = 0
        params = dict(
            muSc=self.muSc,
            mu=self.mu,
            delta=self.delta,
            alpha=self.alpha,
            t=self.t,
            barrier=self.barrier,
            M=self.M,
            added_sinusoid=self.added_sinusoid,
            stagger_ratio=self.stagger_ratio,
            barrier_length=self.barrier_length,
            gfactor=self.gfactor,
            bohr_magneton=self.bohr_magneton,
            hopping_distance=self.hopping_distance,
        )
        for i in tqdm(range(np.size(energies)), desc="Cond",):
            cond = []
            energy = energies[i]
            for b in bValues:
                params["B"] = b
                smatrix = kwant.smatrix(syst, energy, params=params)
                conduct = (
                    smatrix.submatrix((0, 0), (0, 0)).shape[0]  # N
                    - smatrix.transmission((0, 0), (0, 0))  # R_ee
                    + smatrix.transmission(
                        (0, 1), (0, 0)
                    )  # maybe I need to update the 1 to a lattice_constant?
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

    def plot(self):
        syst = NISIN(
            width=self.width,
            noMagnets=self.noMagnets,
            barrier_length=self.barrier_length,
            hopping_distance=self.hopping_distance,
        )

        length = 8 * self.noMagnets - 2 + 2 * self.barrier_length

        return kwant.plotter.plot(
            syst,
            show=False,
            unit="nn",
            site_size=0.20,
            site_color=lambda s: "y"
            if barrier_region(s, self.barrier_length, length, self.width)
            else "b",
        )


def main():
    pass


if __name__ == "__main__":
    main()
