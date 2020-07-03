import kwant
import numpy as np
import tinyarray as ta

s0 = np.identity(2)
sZ = np.array([[1.0, 0.0], [0.0, -1.0]])
sX = np.array([[0.0, 1.0], [1.0, 0.0]])
sY = np.array([[0.0, -1j], [1j, 0.0]])

tauZsig0 = ta.array(np.kron(sZ, s0))
tauXsig0 = ta.array(np.kron(sX, s0))
tauYsig0 = ta.array(np.kron(sY, s0))
tau0sigZ = ta.array(np.kron(s0, sZ))
tau0sigX = ta.array(np.kron(s0, sX))
tau0sigY = ta.array(np.kron(s0, sY))
tauZsigX = ta.array(np.kron(sZ, sX))
tauZsigY = ta.array(np.kron(sZ, sY))
tauYsigY = ta.array(np.kron(sY, sY))
tauYsigZ = ta.array(np.kron(sY, sZ))


def magnetic_phase(position, p):
    counter = np.mod(position - (1 + p["barrier_length"]) * p["hopping_distance"], p["period"])
    theta = (counter / p["period"]) * 2 * np.pi
    return theta


# No need for two t terms, Builder assumes hermiticity
def hopX(site0, site1, p):
    return -p["t"] * tauZsig0 + 1j * p["alpha"] * tauZsigY


def hopY(site0, site1, p):
    return -p["t"] * tauZsig0 - 1j * p["alpha"] * tauZsigX


def sinuB(position, p, stagger_ratio):
    theta = magnetic_phase(position, p)
    return tau0sigY * np.cos(theta) + tau0sigX * np.sin(theta)


def energy_nanomagnet(position, p):
    return 0.5 * p["gfactor"] * p["bohr_magneton"] * p["M"] * sinuB(position, p, p["stagger_ratio"])


# Particle-hole symmetry (lead=True)             {tauY,tauX} {sigY,sigZ}
# Particle-hole symmetry (lead=False & sinusoid) {tauY,tauX} {sigZ}
def onsiteNormal(site, p, lead=False):
    H_no_magnets = (4 * p["t"] - p["mu"]) * tauZsig0 + 0.5 * p["gfactor"] * p["bohr_magneton"] * p["B"] * tau0sigX
    if p["added_sinusoid"] and not lead:
        return H_no_magnets + 0.5 * p["gfactor"] * p["bohr_magneton"] * p["M"] * sinuB(
            site.pos[0], p, p["stagger_ratio"]
        )
    else:
        return H_no_magnets


def onsite_normal(site, p):
    return onsiteNormal(site, p, lead=False)


def onsiteSc(site, p):
    return onsiteNormal(site, p, lead=False) + p["delta"] * tauXsig0


def onsite_lead(site, p):
    return onsiteNormal(site, p, lead=True)


def onsite_barrier(site, p):
    return (
        4 * p["t"]
        - p["mu"]
        + barrier_height_func(p["barrier_height"], p["barrier_length"], p["wire_length"], p["wire_width"], site)
    ) * tauZsig0


def barrier_height_func(barrier_height, barrier_length, wire_length, wire_width, site):
    if barrier_region(site, barrier_length, wire_length, wire_width):
        i = site[1][1] // wire_width
        distance_from_centre = np.abs(((wire_length - 1) / 2) - i)
        dbarrier_height = barrier_height / barrier_length
        negative_distance_from_end = distance_from_centre - ((wire_length - 1) / 2)
        height = barrier_height + dbarrier_height * negative_distance_from_end
        return height
    else:
        raise IndexError("barrier_height called outside of barrier region")


def barrier_region(site, barrier_length, wire_length, wire_width):
    if type(site) == kwant.builder.Site:
        i = site[1][0]
        j = site[1][1]
    else:
        i = site // wire_width
        j = site % wire_width
    return ((0 <= i < barrier_length) or (wire_length - barrier_length <= i < wire_length)) and (0 <= j < wire_width)


def make_system(
    p, onsite_wire=onsiteSc, onsite_barrier=onsite_barrier, onsite_lead=onsite_lead, hopX=hopX, hopY=hopY, norbs=4
):
    def make_wire_and_barriers(p=p, norbs=norbs, onsite_wire=onsite_wire, onsite_barrier=onsite_barrier):
        syst = kwant.Builder()
        lat = kwant.lattice.square(a=p["hopping_distance"], norbs=norbs)

        syst[
            (
                lat(i, j)
                for i in range(p["barrier_length"], p["wire_length"] - p["barrier_length"])
                for j in range(p["wire_width"])
            )
        ] = onsite_wire
        # fmt: off
        syst[
            (
                lat(i, j)
                for i in range(0, p["barrier_length"])
                for j in range(p["wire_width"])
            )
        ] = onsite_barrier
        # fmt: on
        syst[
            (
                lat(i, j)
                for i in range(p["wire_length"] - p["barrier_length"], p["wire_length"])
                for j in range(p["wire_width"])
            )
        ] = onsite_barrier

        # Hopping:
        syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
        syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY

        return syst

    def make_lead(
        p=p, onsite_lead=onsite_lead, lead_hopX=-p["t"] * tauZsig0, lead_hopY=-p["t"] * tauZsig0, norbs=norbs
    ):
        # Conservation law must separate the electron-hole degree of freedom -> tauZ
        # Particle-hole symmetry operator must be figured out for the Hamiltonian
        lead = kwant.Builder(
            kwant.TranslationalSymmetry((-p["hopping_distance"], 0)), conservation_law=-tauZsig0, particle_hole=tauYsigZ
        )
        lat = kwant.lattice.square(a=p["hopping_distance"], norbs=norbs)
        lead[(lat(0, j) for j in range(p["wire_width"]))] = onsite_lead
        lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = lead_hopX
        lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = lead_hopY
        return lead

    syst = make_wire_and_barriers()
    lead = make_lead()
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def NISIN(params):
    return make_system(params).finalized()


def NININ(params):
    return make_system(params, onsite_wire=onsite_normal).finalized()
