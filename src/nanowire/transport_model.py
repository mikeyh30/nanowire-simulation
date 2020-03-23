import kwant
import numpy as np
import tinyarray as ta

s0 = np.identity(2)
sZ = np.array([[1.0, 0.0], [0.0, -1.0]])
sX = np.array([[0.0, 1.0], [1.0, 0.0]])
sY = np.array([[0.0, -1j], [1j, 0.0]])

tauZ = ta.array(np.kron(sZ, s0))
tauX = ta.array(np.kron(sX, s0))
tauY = ta.array(np.kron(sY, s0))
sigZ = ta.array(np.kron(s0, sZ))
sigX = ta.array(np.kron(s0, sX))
sigY = ta.array(np.kron(s0, sY))
tauZsigX = ta.array(np.kron(sZ, sX))
tauZsigY = ta.array(np.kron(sZ, sY))
tauYsigY = ta.array(np.kron(sY, sY))


def magnetic_phase(position, barrier_length, hopping_distance, period):
    counter = np.mod(position - (1 + barrier_length) * hopping_distance, period)
    theta = (counter / period) * 2 * np.pi
    return theta


def hopX(site0, site1, t, alpha):
    return -t * tauZ + 1j * alpha * tauZsigY


def hopY(site0, site1, t, alpha):
    return -t * tauZ - 1j * alpha * tauZsigX


def sinuB(position, barrier_length, hopping_distance, period, stagger_ratio):
    theta = magnetic_phase(position, barrier_length, hopping_distance, period)
    return sigY * np.cos(theta) + sigX * np.sin(theta)


# This is the onsite Hamiltonian, this is where the B-field can be varied.
def onsiteSc(
    site,
    muSc,
    t,
    B,
    delta,
    M,
    added_sinusoid,
    barrier_length,
    stagger_ratio,
    gfactor,
    bohr_magneton,
    period,
    hopping_distance,
):
    if added_sinusoid:
        return (
            (4 * t - muSc) * tauZ
            + 0.5 * gfactor * bohr_magneton * B * sigX
            + delta * tauX
            + 0.5
            * gfactor
            * bohr_magneton
            * M
            * sinuB(
                site.pos[0], barrier_length, hopping_distance, period, stagger_ratio
            )
        )
    else:
        return (
            (4 * t - muSc) * tauZ
            + 0.5 * gfactor * bohr_magneton * B * sigX
            + delta * tauX
        )


def onsiteNormal(site, mu, t):
    return (4 * t - mu) * tauZ


def onsiteBarrier(site, mu, t, barrier):
    return (4 * t - mu + barrier) * tauZ


def make_lead(width, hopping_distance, onsiteH=onsiteNormal, hopX=hopX, hopY=hopY):
    lead = kwant.Builder(
        kwant.TranslationalSymmetry((-hopping_distance, 0)),
        conservation_law=-tauZ,
        particle_hole=tauYsigY,
    )
    lat = kwant.lattice.square(a=hopping_distance, norbs=4)
    lead[(lat(0, j) for j in range(width))] = onsiteH
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    return lead


def barrier_region(site, barrier_length, length, width):
    i = site // width
    j = site % width
    return ((0 <= i < barrier_length) or (length - barrier_length <= i < length)) and (
        0 <= j < width
    )


def make_wire(
    width,
    length,
    barrier_length,
    hopping_distance,
    hamiltonian_wire=onsiteSc,
    hamiltonian_barrier=onsiteBarrier,
    hamiltonian_normal=onsiteNormal,
    hopX=hopX,
    hopY=hopY,
):

    syst = kwant.Builder()
    lat = kwant.lattice.square(a=hopping_distance, norbs=4)

    syst[
        (
            lat(i, j)
            for i in range(barrier_length, length - barrier_length)
            for j in range(width)
        )
    ] = onsiteSc
    # fmt: off
    syst[
        (
            lat(i, j)
            for i in range(0, barrier_length + 1)
            for j in range(width)
        )
    ] = onsiteBarrier
    # fmt: on
    syst[
        (
            lat(i, j)
            for i in range(length - barrier_length, length)
            for j in range(width)
        )
    ] = onsiteBarrier

    # Hopping:
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    lead = make_lead(width, hopping_distance, onsiteNormal, hopX, hopY)
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def NISIN(width=7, length=40, barrier_length=1, hopping_distance=1):
    return make_wire(width, length, barrier_length, hopping_distance).finalized()
