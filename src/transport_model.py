import kwant
import numpy as np
import tinyarray as ta

bohr_magneton =  58E-6 # eVT^-1 # physical_constants['Bohr magneton'][0]
lattice_constant_InAs = 200 # angstroms # 6.0583E-10 # 20E-3 # might need to change this.

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

lat = kwant.lattice.square(a=lattice_constant_InAs,norbs=4)

def hopX(site0, site1, t, alpha):
    return -t * tauZ + 1j * alpha * tauZsigY

def hopY(site0, site1, t, alpha):
    return -t * tauZ - 1j * alpha * tauZsigX

def sinuB(theta, stagger_ratio):
    # ssin, scos = rick_fourier(theta)
    # return sigY*scos + sigX*ssin
    return sigY * np.cos(theta) + sigX * np.sin(theta)
# This is the onsite Hamiltonian, this is where the B-field can be varied.

def onsiteSc(site, muSc, t, B, Delta, M, addedSinu, barrier_length, stagger_ratio):
    gfactor = 10  # should be 10 in the real units
    if addedSinu:
        counter = np.mod(site.pos[0] - 1 - barrier_length, 16)
        if -1 < counter < 4:
            theta = 0
        elif 3 < counter < 8:
            theta = 0.2 * (counter - 3) * np.pi
        elif 7 < counter < 12:
            theta = np.pi
        else:
            theta = 0.2 * (counter - 6) * np.pi

        return (
            (4 * t - muSc) * tauZ
            + 0.5 * gfactor * bohr_magneton * B * sigX
            + Delta * tauX
            + 0.5 * gfactor * bohr_magneton * M * sinuB(theta, stagger_ratio)
        )
    else:
        return (
            (4 * t - muSc) * tauZ
            + 0.5 * gfactor * bohr_magneton * B * sigX
            + Delta * tauX
        )

def onsiteNormal(site, mu, t):
    return (4 * t - mu) * tauZ

def onsiteBarrier(site, mu, t, barrier):
    return (4 * t - mu + barrier) * tauZ

def make_lead(width, onsiteH=onsiteNormal, hopX=hopX, hopY=hopY):
    lead = kwant.Builder(
        kwant.TranslationalSymmetry((-lattice_constant_InAs, 0)),
        conservation_law=-tauZ,
        particle_hole=tauYsigY,
    )
    lat = kwant.lattice.square(a=lattice_constant_InAs, norbs=4)
    lead[(lat(0, j) for j in range(width))] = onsiteH
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopX
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopY
    return lead

def barrier_region(site, barrier_length, length, width):
    i = site//width
    j = site%width
    return(
        (
            (0 <= i < barrier_length)
            or
            (length - barrier_length <= i < length)
        )
        and
        (
            0 <= j < width
        )
      )

def make_wire(width, length, barrier_length, hamiltonian_wire=onsiteSc,
        hamiltonian_barrier=onsiteBarrier, hamiltonian_normal=onsiteNormal,
        hopX=hopX, hopY=hopY):

    syst=kwant.Builder()

    syst[
        (
            lat(i, j)
            for i in range(barrier_length, length - barrier_length)
            for j in range(width)
        )
    ] = onsiteSc

    syst[
        (
            lat(i, j)
            for i in range(0, barrier_length + 1)
            for j in range(width)
        )
    ] = onsiteBarrier

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
    lead = make_lead(width, onsiteNormal, hopX, hopY)
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def NISIN(width=7, noMagnets=5, barrier_length=1):

    length = 8 * noMagnets - 2 + 2 * barrier_length  # "2*(noMagnets - 1)"
    return make_wire(width, length, barrier_length).finalized()