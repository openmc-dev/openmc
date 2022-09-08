from math import pi

import openmc
import openmc.stats
import numpy as np
from pytest import approx


def test_source():
    space = openmc.stats.Point()
    energy = openmc.stats.Discrete([1.0e6], [1.0])
    angle = openmc.stats.Isotropic()

    src = openmc.Source(space=space, angle=angle, energy=energy)
    assert src.space == space
    assert src.angle == angle
    assert src.energy == energy

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert elem.find('space') is not None
    assert elem.find('angle') is not None
    assert elem.find('energy') is not None

    src = openmc.Source.from_xml_element(elem)
    assert isinstance(src.angle, openmc.stats.Isotropic)
    assert src.space.xyz == [0.0, 0.0, 0.0]
    assert src.energy.x == [1.0e6]
    assert src.energy.p == [1.0]
    assert src.strength == 1.0


def test_spherical_uniform():
    r_outer = 2.0
    r_inner = 1.0
    thetas = (0.0, pi/2)
    phis = (0.0, pi)
    origin = (0.0, 1.0, 2.0)

    sph_indep_function = openmc.stats.spherical_uniform(r_outer,
                                                        r_inner,
                                                        thetas,
                                                        phis,
                                                        origin)

    assert isinstance(sph_indep_function, openmc.stats.SphericalIndependent)


def test_source_file():
    filename = 'source.h5'
    src = openmc.Source(filename=filename)
    assert src.file == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib


def test_source_dlopen():
    library = './libsource.so'
    src = openmc.Source(library=library)
    assert src.library == library

    elem = src.to_xml_element()
    assert 'library' in elem.attrib


def test_source_xml_roundtrip():
    # Create a source and write to an XML element
    space = openmc.stats.Box([-5., -5., -5.], [5., 5., 5.])
    energy = openmc.stats.Discrete([1.0e6, 2.0e6, 5.0e6], [0.3, 0.5, 0.2])
    angle = openmc.stats.PolarAzimuthal(
        mu=openmc.stats.Uniform(0., 1.),
        phi=openmc.stats.Uniform(0., 2*pi),
        reference_uvw=(0., 1., 0.)
    )
    src = openmc.Source(
        space=space, angle=angle, energy=energy,
        particle='photon', strength=100.0
    )
    elem = src.to_xml_element()

    # Read from XML element and make sure data is preserved
    new_src = openmc.Source.from_xml_element(elem)
    assert isinstance(new_src.space, openmc.stats.Box)
    np.testing.assert_allclose(new_src.space.lower_left, src.space.lower_left)
    np.testing.assert_allclose(new_src.space.upper_right, src.space.upper_right)
    assert isinstance(new_src.energy, openmc.stats.Discrete)
    np.testing.assert_allclose(new_src.energy.x, src.energy.x)
    np.testing.assert_allclose(new_src.energy.p, src.energy.p)
    assert isinstance(new_src.angle, openmc.stats.PolarAzimuthal)
    assert new_src.angle.mu.a == src.angle.mu.a
    assert new_src.angle.mu.b == src.angle.mu.b
    assert new_src.angle.phi.a == src.angle.phi.a
    assert new_src.angle.phi.b == src.angle.phi.b
    np.testing.assert_allclose(new_src.angle.reference_uvw, src.angle.reference_uvw)
    assert new_src.particle == src.particle
    assert new_src.strength == approx(src.strength)
