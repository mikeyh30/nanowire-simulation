from pytest import approx, raises
from nanowire.transport_model import * 
import kwant
import tinyarray as ta

p = {
    "wire_width": 7,
    "wire_length": 80,
    "barrier_length": 1,
    "stagger_ratio": 0.5,
    "period": 16,
    "M": 1,
    "m_max": 5,
    "hopping_distance": 1,
    "added_sinusoid": True,
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


def site(x,y):
    fam = kwant.lattice.general(ta.identity(2))
    return kwant.builder.Site(fam,(x,y))


def test_magnetic_phase():
    assert magnetic_phase(2,p) == 0
    assert magnetic_phase(0,p) == approx(5.4977,0.0001)
    with raises(TypeError):
        magnetic_phase("a",p)


def test_barrier_region():
    assert barrier_region(p,site=site(0,1))
    assert not barrier_region(p,site=site(2,4))
    assert barrier_region(p,x=0,y=4)
    assert not barrier_region(p,x=0,y=9)
    with raises(TypeError):
        barrier_region(p,site=site("s",1))
    with raises(TypeError):
        barrier_region(p,x=(0,1),y=2)
    with raises(TypeError):
        barrier_region(p,site)
    with raises(TypeError):
        barrier_region(p,ste=site("s",1))


def test_wire_region():
    assert wire_region(p,site=site(0,1))
    assert wire_region(p,site=site(2,4))
    assert wire_region(p,x=0,y=4)
    assert not wire_region(p,x=0,y=9)
    assert not wire_region(p,x=-2,y=2)
    with raises(TypeError):
        wire_region(p,site=site("s",1))
    with raises(TypeError):
        wire_region(p,x=(0,1),y=2)
    with raises(TypeError):
        wire_region(p,site)
    with raises(TypeError):
        wire_region(p,ste=site("s",1))


def test_sc_region():
    assert not sc_region(p,site=site(0,1))
    assert sc_region(p,site=site(3,5))
    assert not sc_region(p,x=0,y=4)
    assert not sc_region(p,x=0,y=9)
    assert not sc_region(p,x=80,y=4)
    with raises(TypeError):
        sc_region(p,site=site("s",1))
    with raises(TypeError):
        sc_region(p,x=(0,1),y=2)
    with raises(TypeError):
        sc_region(p,site)
    with raises(TypeError):
        sc_region(p,ste=site("s",1))


def test_f_M():
    assert f_M(2,p) == (0,1)
    assert f_M(0,p) == approx((-0.70711,0.70711),0.0001)
    with raises(TypeError):
        f_M([0,1])
    



