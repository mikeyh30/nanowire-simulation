from tqdm import tqdm
import kwant
import tinyarray as ta
import numpy as np
import scipy.sparse.linalg
from nanowire.nanomagnet_field import rick_fourier
from nanowire.transport_model import NISIN, barrier_region, magnetic_phase
import matplotlib.pyplot as plt


class Nanowire:
    def __init__(
        self,
        parameters
    ):
        parameters['t'] = 3.83 / (parameters['effective_mass'] * (parameters['hopping_distance'] ** 2))
        parameters['alpha'] = parameters['alpha_R'] / parameters['hopping_distance']
        self.parameters = parameters

    def spectrum(self, bValues=np.linspace(0, 1.0, 201)):
        syst = NISIN(self.parameters)
        energies = []
        critB = 0
        for b in tqdm(bValues, desc="Spec",):
            self.parameters["B"] = b
            H = syst.hamiltonian_submatrix(sparse=True, params=self.parameters)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)
            if critB == 0 and np.abs(eigs[10] - eigs[9]) / 2 < 1e-4:
                critB = b

        outcome = dict(B=bValues, E=energies, CritB=critB)
        return outcome
    
    def paper_spectrum(self, m_values):
        syst = NISIN(self.parameters)
        energies = []
        critM = 0
        for m in tqdm(m_values, desc="Spec",):
            self.parameters["M"] = m
            H = syst.hamiltonian_submatrix(sparse=True, params=self.parameters)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)
            if critM == 0 and np.abs(eigs[10] - eigs[9]) / 2 < 1e-4:
                critM = m

        outcome = dict(M=m_values, E=energies, CritM=critM)
        return outcome
    
    def conductances(
        self,
        bValues=np.linspace(0, 1.0, 201),
        energies=[1e-6 * i for i in range(-120, 120)],
    ):
        syst = NISIN(self.parameters)
        data = []
        critB = 0
        for energy in tqdm(energies, desc="Cond",):
            cond = []
            for b in bValues:
                self.parameters["B"] = b
                smatrix = kwant.smatrix(syst, energy, params=self.parameters)
                conduct = (
                    smatrix.submatrix((0, 0), (0, 0)).shape[0]  # N
                    - smatrix.transmission((0, 0), (0, 0))  # R_ee
                    + smatrix.transmission((0, 1), (0, 0))  # R_he
                )
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

    def plot(self, ax_model, ax_x, ax_y):
        syst = NISIN(self.parameters)
        length_A = syst.pos(self.parameters['wire_width'] * self.parameters['wire_length'] - 1)[0]
        array_A = np.arange(length_A)
        phi = magnetic_phase(
            array_A, self.parameters['barrier_length'],
            self.parameters['hopping_distance'], self.parameters['period']
        )
        ax_x.plot(array_A, self.parameters['M'] * np.sin(phi))
        ax_y.plot(array_A, self.parameters['M'] * np.cos(phi))

        return kwant.plotter.plot(
            syst,
            show=False,
            unit="nn",
            site_size=0.20,
            site_color=lambda s: "y"
            if barrier_region(s, self.parameters['barrier_length'],
                              self.parameters['wire_length'],
                              self.parameters['wire_width'])
            else "b",
            ax=ax_model,
        )


def main():
    pass


if __name__ == "__main__":
    main()
