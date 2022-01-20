import warnings
import kwant
warnings.filterwarnings("ignore", category=DeprecationWarning)
import kwant.continuum
warnings.filterwarnings("ignore", category=DeprecationWarning)
import numpy as np
import tinyarray as ta

s0 = np.identity(2)
sZ = np.array([[1.0, 0.0], [0.0, -1.0]])
sX = np.array([[0.0, 1.0], [1.0, 0.0]])
sY = np.array([[0.0, -1j], [1j, 0.0]])


sigma = {"0": s0, "X": sX, "Y": sY, "Z": sZ}


def hamiltonian(rashba=False, superconducting=False, micromagnets=False):
    hamiltonian_normal = (
        "A * (k_x**2 + k_y**2) * kron(sigma_z, sigma_0) +"
        "-mu_wire * kron(sigma_z, sigma_0) + "
        "0.5 * gfactor * bohr_magneton * B * kron(sigma_0, sigma_x)"
    )

    hamiltonian_rashba = "+ alpha_R * k_x * kron(sigma_z, sigma_y) + alpha_R * k_y * kron(sigma_z, sigma_x)"

    hamiltonian_superconducting = "+ delta * kron(sigma_x, sigma_0)"

    hamiltonian_micromagnets_x = "+ 0.5 * gfactor * bohr_magneton * M * sinM(x) * kron(sigma_0, sigma_x)"

    hamiltonian_micromagnets_y = "+ 0.5 * gfactor * bohr_magneton * M * cosM(y) * kron(sigma_0, sigma_y)"

    hamiltonian = (
        hamiltonian_normal
        + rashba * hamiltonian_rashba
        + superconducting * hamiltonian_superconducting
        + micromagnets * hamiltonian_micromagnets_x
        + micromagnets * hamiltonian_micromagnets_y
    )

    return hamiltonian


def norbitals(p):
    if p["delta"] == 0:
        return 2, False
    elif p["delta"] > 0:
        return 4, True
    else:
        raise ValueError


def magnetic_phase(position, p):
    counter = np.mod(position - (1 + p["barrier_length"]) * p["hopping_distance"], p["period"])
    theta = (counter / p["period"]) * 2 * np.pi
    return theta

    
def barrier_height_func(barrier_height, barrier_length, wire_length, wire_width, site):
    if barrier_region(site, wire_width, wire_length, barrier_length):
        i = site[1][1] // wire_width
        distance_from_centre = np.abs(((wire_length - 1) / 2) - i)
        dbarrier_height = barrier_height / barrier_length
        negative_distance_from_end = distance_from_centre - ((wire_length - 1) / 2)
        height = barrier_height + dbarrier_height * negative_distance_from_end
        return height
    else:
        raise IndexError("barrier_height called outside of barrier region")


def central_region(site, w, c_l, b_l, a):
    (x, y) = site.pos
    return b_l*a <= x < c_l*a + b_l*a and 0 <= y < w*a


def barrier_region(site, w, c_l, b_l, a):
    if type(site) == kwant.builder.Site:
        (x, y) = site.pos
    else:
        x = site // w
        y = site % w
    left_barrier = 0 <= x < b_l*a and 0 <= y < w*a
    right_barrier = c_l*a + b_l*a <= x < c_l*a + 2 * b_l*a and 0 <= y < w*a
    return left_barrier or right_barrier


def make_system(
#    p, onsite_wire=onsite_wire_sc, onsite_barrier=onsite_barrier, onsite_lead=onsite_lead, hopX=hopX, hopY=hopY, norbs=4
    w, c_l, b_l, central_ham, barrier_ham, a
):
    syst = kwant.Builder()

    wrapped_central_region = lambda site : central_region(site, w, c_l, b_l, a)
    wrapped_barrier_region = lambda site : barrier_region(site, w, c_l, b_l, a)

    syst.fill(barrier_ham, wrapped_barrier_region, (0, 0))
    syst.fill(central_ham, wrapped_central_region, (b_l*a, 0))
    syst.fill(barrier_ham, wrapped_barrier_region, (b_l*a + c_l*a, 0))

    def make_lead():#p=p, onsite_lead=onsite_lead, lead_hopX=lead_hopX, lead_hopY=lead_hopY, norbs=norbs):
        # Conservation law must separate the electron-hole degree of freedom -> tauZ
        # Particle-hole symmetry operator must be figured out for the Hamiltonian
        tau_z = np.kron(np.array([[1, 0], [0, -1]]), np.eye(2))
        tau_y_sig_z = np.kron(np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]]))

        lead = kwant.Builder(
            kwant.TranslationalSymmetry((-a, 0)),
            conservation_law=-tau_z,
            particle_hole=tau_y_sig_z,
        )

        def lead_shape(site):
            (x, y) = site.pos
            return 0 <= x < 1 and 0 <= y < w*a

        lead.fill(barrier_ham, lead_shape, (0, 0))
        return lead

    # TODO Look into where the different chemical potentials are used for each lead
    l_lead = make_lead()#onsite_lead=onsite_l_lead)
    r_lead = make_lead()#onsite_lead=onsite_r_lead)
    syst.attach_lead(l_lead)
    syst.attach_lead(r_lead.reversed())

    return syst


def NIXIN(params):
    norbs, superconducting = norbitals(params)
    w = params['wire_width']
    c_l = params['wire_length'] - 2 * params['barrier_length']
    b_l = params['barrier_length']
    a = params['hopping_distance']
    central_ham = kwant.continuum.discretize(hamiltonian(True, superconducting, True), coords="xy", grid=a)
    barrier_ham = kwant.continuum.discretize(hamiltonian(True, False, True), coords="xy", grid=a)
    return make_system(w, c_l, b_l, central_ham, barrier_ham, a).finalized()
