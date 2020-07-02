from pytest import approx
from nanowire.nanowire import *
import math

p = {
    "wire_width": 7,
    "wire_length": 80,
    "barrier_length": 1,
    "stagger_ratio": 0.5,
    "period": 16,
    "M": 0,
    "m_max": 5,
    "hopping_distance": 1,
    "added_sinusoid": False,
    "B": 1,
    "b_max": 3,
    "bohr_magneton": 1,
    "alpha_R": 0.32,
    "delta": 0.1,
    "gfactor": 1,
    "effective_mass": 1,
    "muSc": 0.22,
    "mu": 0.3,
    "barrier_height": 2,
}


def test_spectrum():
    nanowire = Nanowire(p)
    outcome = nanowire.spectrum(B_values=np.linspace(0, 10, 3), tolerance=0.001)
    assert outcome["CritB"] == 5.0
    outcome2 = nanowire.spectrum(B_values=np.linspace(0, 10, 3), tolerance=0.00001)
    assert math.isnan(outcome2["CritB"])
    return outcome


def test_spin_density():
    nanowire = Nanowire(p)
    spin_z, _ = nanowire.spin_density(0, 1, True)
    assert spin_z == approx(0.4766, 0.0001)


if __name__ == "__main__":
    test_spin_density()
