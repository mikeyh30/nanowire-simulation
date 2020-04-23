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


def magnetic_phase(position, p):
    counter = np.mod(position - (1 + p['barrier_length']) * p['hopping_distance'], p['period'])
    theta = (counter / p['period']) * 2 * np.pi
    return theta


# No need for two t terms, Builder assumes hermiticity
def hopX(site0, site1, p):
    return -p['t'] * tauZsig0 + 1j * p['alpha'] * tauZsigY


def hopY(site0, site1, p):
    return -p['t'] * tauZsig0 - 1j * p['alpha'] * tauZsigX


def sinuB(position, p, stagger_ratio):
    theta = magnetic_phase(position, p)
    return tau0sigY * np.cos(theta) + tau0sigX * np.sin(theta)


def energy_nanomagnet(position,p):
    return (
        0.5
        * p['gfactor']
        * p['bohr_magneton']
        * p['M']
        * sinuB(position, p, p['stagger_ratio'])
    )


# This is the onsite Hamiltonian, this is where the B-field can be varied.
def onsiteSc(site,p):
    H_no_magnets = (
        (4 * p['t'] - p['mu']) * tauZsig0
        + 0.5 * p['gfactor'] * p['bohr_magneton'] * p['B'] * tau0sigX
        + p['delta'] * tauXsig0
    )
    if p['added_sinusoid']:  # Might consider changing this to if M:, if float zero is good
        return (
            H_no_magnets
            + 0.5 * p['gfactor'] * p['bohr_magneton'] * p['M'] * sinuB(site.pos[0], p, p['stagger_ratio'])
        )
    else:
        return H_no_magnets


def onsiteNormal(site, p):
    return (4 * p['t'] - p['mu']) * tauZsig0


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


def onsiteBarrier(site, p):
    return (4 * p['t'] - p['mu'] + \
        barrier_height_func(p['barrier_height'], p['barrier_length'], p['wire_length'], p['wire_width'], site)) * tauZsig0


def make_lead(p, onsiteH=onsiteNormal, hopX=hopX, hopY=hopY):
    lead = kwant.Builder(
        kwant.TranslationalSymmetry((-p['hopping_distance'], 0)),
        conservation_law=-tauZsig0,
        particle_hole=tauYsigY,
    )
    lat = kwant.lattice.square(a=p['hopping_distance'], norbs=4)
    lead[(lat(0, j) for j in range(p['wire_width']))] = onsiteH
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = -p['t'] * tauZsig0
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = -p['t'] * tauZsig0
    return lead


def barrier_region(site, barrier_length, wire_length, wire_width):
    if type(site) == kwant.builder.Site:
        i = site[1][0]
        j = site[1][1]
    else:
        i = site // wire_width
        j = site % wire_width
    return (
        (0 <= i < barrier_length) or (wire_length - barrier_length <= i < wire_length)
    ) and (0 <= j < wire_width)


def make_wire(
    params,
    hamiltonian_wire=onsiteSc,
    hamiltonian_barrier=onsiteBarrier,
    hamiltonian_normal=onsiteNormal,
    hopX=hopX,
    hopY=hopY,
):
    wire_width = params["wire_width"]
    wire_length = params["wire_length"]
    barrier_length = params["barrier_length"]
    hopping_distance = params["hopping_distance"]

    syst = kwant.Builder()
    lat = kwant.lattice.square(a=hopping_distance, norbs=4)

    syst[
        (
            lat(i, j)
            for i in range(barrier_length, wire_length - barrier_length)
            for j in range(wire_width)
        )
    ] = onsiteSc
    # fmt: off
    syst[
        (
            lat(i, j)
            for i in range(0, barrier_length)
            for j in range(wire_width)
        )
    ] = onsiteBarrier
    # fmt: on
    syst[
        (
            lat(i, j)
            for i in range(wire_length - barrier_length, wire_length)
            for j in range(wire_width)
        )
    ] = onsiteBarrier

    # Hopping:
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    lead = make_lead(params, onsiteNormal, hopX, hopY)
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def NISIN(params):
    return make_wire(params).finalized()
