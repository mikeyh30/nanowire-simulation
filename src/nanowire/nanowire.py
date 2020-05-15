from tqdm import tqdm
import kwant
import numpy as np
import scipy.sparse.linalg
from nanowire.transport_model import NISIN, barrier_region, magnetic_phase


def find_critical_field(
    B_values, energies, min_topological_gap, tolerance=1e-5, energy_level_count=20
):
    topological_B = []
    topological_gap = []
    middle_energy = energy_level_count // 2
    for bidx, b in enumerate(B_values):
        if np.abs(energies[bidx][middle_energy] - energies[bidx][middle_energy - 1]) / 2 < tolerance:
            if np.abs(energies[bidx][middle_energy + 1] - energies[bidx][middle_energy]) > min_topological_gap:
                topological_B.append(b)
                topological_gap.append(np.abs(energies[bidx][middle_energy + 1]- energies[bidx][middle_energy]))
    return topological_B, topological_gap


class Nanowire:
    def __init__(self, parameters):
        parameters["t"] = 3.83 / (parameters["effective_mass"] * (parameters["hopping_distance"] ** 2))
        parameters["alpha"] = parameters["alpha_R"] / parameters["hopping_distance"]
        self.parameters = parameters

    def spectrum(self, simulation_run, B_values=np.linspace(0, 1.0, 201)):
        syst = NISIN(self.parameters)
        energies = []
        critB = 0
        for b in tqdm(B_values, desc=simulation_run + "/spec",):
            self.parameters["B"] = b
            newparams = {}
            newparams['p'] = self.parameters
            H = syst.hamiltonian_submatrix(sparse=True, params=newparams)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)

        topological_B_values, topological_gap = find_critical_field(B_values, energies, 0.15 * self.parameters["delta"])
        if not topological_B_values:
            topological_B_values.append(np.nan)

        outcome = dict(
            B=B_values,
            E=energies,
            CritB=topological_B_values[0],
            topological_B_values=topological_B_values,
            topological_gap=topological_gap,
        )
        return outcome

    def magnetization_spectrum(self, simulation_run, M_values):
        syst = NISIN(self.parameters)
        energies = []
        critM = 0
        for m in tqdm(M_values, desc=simulation_run + "/mag-spec",):
            self.parameters["M"] = m
            H = syst.hamiltonian_submatrix(sparse=True, params=self.parameters)
            H = H.tocsc()
            # k is the number of eigenvalues, and find them near sigma.
            eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
            eigs = np.sort(eigs[0])
            energies.append(eigs)

        topological_M_values = find_critical_field(M_values, energies, 0.15 * self.parameters["delta"])

        outcome = dict(M=M_values, E=energies, CritM=topological_M_values[0])
        return outcome

    def conductances(
        self,
        simulation_run,
        B_values=np.linspace(0, 1.0, 201),
        energies=[1e-6 * i for i in range(-120, 120)],
    ):
        syst = NISIN(self.parameters)
        data = []
        critB = 0
        for energy in tqdm(energies, desc=simulation_run + "/cond",):
            cond = []
            for b in B_values:
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

        outcome = dict(B=B_values, BiasV=energies, Cond=data, CritB=critB)
        return outcome

    def plot(self, ax_model, ax_x, ax_y):
        syst = NISIN(self.parameters)
        length_A = syst.pos(self.parameters["wire_width"] * self.parameters["wire_length"] - 1)[0]
        array_A = np.arange(length_A)
        phi = magnetic_phase(array_A, self.parameters)
        ax_x.plot(array_A, self.parameters["M"] * np.sin(phi))
        ax_y.plot(array_A, self.parameters["M"] * np.cos(phi))

        return kwant.plotter.plot(
            syst,
            show=False,
            unit="nn",
            site_size=0.20,
            site_color=lambda s: "y"
            if barrier_region(
                s,
                self.parameters["barrier_length"],
                self.parameters["wire_length"],
                self.parameters["wire_width"],
            )
            else "b",
            ax=ax_model,
        )


def main():
    pass


if __name__ == "__main__":
    main()
