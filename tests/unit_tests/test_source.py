from math import pi

import openmc
import openmc.lib
import openmc.stats
import numpy as np
import pytest
from pytest import approx

from tests.regression_tests import config


def test_source():
    space = openmc.stats.Point()
    energy = openmc.stats.Discrete([1.0e6], [1.0])
    angle = openmc.stats.Isotropic()

    src = openmc.IndependentSource(space=space, angle=angle, energy=energy)
    assert src.space == space
    assert src.angle == angle
    assert src.energy == energy

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert elem.find('space') is not None
    assert elem.find('angle') is not None
    assert elem.find('energy') is not None

    src = openmc.IndependentSource.from_xml_element(elem)
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

def test_point_cloud():
    point_list = [[1,0,0], [0,1,0], [0,0,1]]
    positions = np.asarray(point_list)
    strengths = [1,2,3]

    space = openmc.stats.PointCloud(positions, strengths)
    np.testing.assert_equal(space.positions, positions)
    np.testing.assert_equal(space.strengths, strengths)

    space = openmc.stats.PointCloud(point_list, strengths)
    np.testing.assert_equal(space.positions, positions)
    np.testing.assert_equal(space.strengths, strengths)

    energy = openmc.stats.Discrete([1.0e6], [1.0])
    angle = openmc.stats.Isotropic()

    src = openmc.IndependentSource(space=space, angle=angle, energy=energy)
    assert src.space == space
    np.testing.assert_equal(src.space.positions, positions)
    np.testing.assert_equal(src.space.strengths, strengths)

    elem = src.to_xml_element()
    src = openmc.IndependentSource.from_xml_element(elem)
    np.testing.assert_equal(src.space.positions, positions)
    np.testing.assert_equal(src.space.strengths, strengths)


def test_source_file():
    filename = 'source.h5'
    src = openmc.FileSource(path=filename)
    assert src.path == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib


def test_source_dlopen():
    library = './libsource.so'
    src = openmc.CompiledSource(library=library)
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
    src = openmc.IndependentSource(
        space=space, angle=angle, energy=energy,
        particle='photon', strength=100.0
    )
    elem = src.to_xml_element()

    # Read from XML element and make sure data is preserved
    new_src = openmc.IndependentSource.from_xml_element(elem)
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


@pytest.fixture
def sphere_box_model():
    # Model with two spheres inside a box
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    sph1 = openmc.Sphere(x0=3, r=1.0)
    sph2 = openmc.Sphere(x0=-3, r=1.0)
    cube = openmc.model.RectangularParallelepiped(
        -5., 5., -5., 5., -5., 5., boundary_type='reflective'
    )
    cell1 = openmc.Cell(fill=mat, region=-sph1)
    cell2 = openmc.Cell(fill=mat, region=-sph2)
    non_source_region = +sph1 & +sph2 & -cube
    cell3 = openmc.Cell(region=non_source_region)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2, cell3])
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'

    return model, cell1, cell2, cell3


def test_constraints_independent(sphere_box_model, run_in_tmpdir):
    model, cell1, cell2, cell3 = sphere_box_model

    # Set up a box source with rejection on the spherical cell
    space = openmc.stats.Box((-4., -1., -1.), (4., 1., 1.))
    model.settings.source = openmc.IndependentSource(
        space=space, constraints={'domains': [cell1, cell2]}
    )

    # Load up model via openmc.lib and sample source
    model.export_to_model_xml()
    openmc.lib.init()
    particles = openmc.lib.sample_external_source(1000)

    # Make sure that all sampled sources are within one of the spheres
    for p in particles:
        assert p.r in (cell1.region | cell2.region)
        assert p.r not in cell3.region

    openmc.lib.finalize()


def test_constraints_mesh(sphere_box_model, run_in_tmpdir):
    model, cell1, cell2, cell3 = sphere_box_model

    bbox = cell3.bounding_box
    mesh = openmc.RegularMesh()
    mesh.lower_left = bbox.lower_left
    mesh.upper_right = bbox.upper_right
    mesh.dimension = (2, 1, 1)

    left_source = openmc.IndependentSource()
    right_source = openmc.IndependentSource()
    model.settings.source = openmc.MeshSource(
        mesh, [left_source, right_source], constraints={'domains': [cell1, cell2]}
    )

    # Load up model via openmc.lib and sample source
    model.export_to_model_xml()
    openmc.lib.init()
    particles = openmc.lib.sample_external_source(1000)

    # Make sure that all sampled sources are within one of the spheres
    for p in particles:
        assert p.r in (cell1.region | cell2.region)
        assert p.r not in cell3.region

    openmc.lib.finalize()


def test_constraints_file(sphere_box_model, run_in_tmpdir):
    model = sphere_box_model[0]

    # Create source file with randomly sampled source sites
    rng = np.random.default_rng()
    energy = rng.uniform(0., 1e6, 10_000)
    time = rng.uniform(0., 1., 10_000)
    particles = [openmc.SourceParticle(E=e, time=t) for e, t in zip(energy, time)]
    openmc.write_source_file(particles, 'uniform_source.h5')

    # Use source file
    model.settings.source = openmc.FileSource(
        'uniform_source.h5',
        constraints={
            'time_bounds': [0.25, 0.75],
            'energy_bounds': [500.e3, 1.0e6],
        }
    )

    # Load up model via openmc.lib and sample source
    model.export_to_model_xml()
    openmc.lib.init()
    particles = openmc.lib.sample_external_source(1000)

    # Make sure that all sampled sources are within energy/time bounds
    for p in particles:
        assert 0.25 <= p.time <= 0.75
        assert 500.e3 <= p.E <= 1.0e6

    openmc.lib.finalize()


@pytest.mark.skipif(config['mpi'], reason='Not compatible with MPI')
def test_rejection_limit(sphere_box_model, run_in_tmpdir):
    model, cell1 = sphere_box_model[:2]

    # Define a point source that will get rejected 100% of the time
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((-3., 0., 0.)),
        constraints={'domains': [cell1]}
    )

    # Confirm that OpenMC doesn't run in an infinite loop. Note that this may
    # work when running with MPI since it won't necessarily capture the error
    # message correctly
    with pytest.raises(RuntimeError, match="rejected"):
        model.run(openmc_exec=config['exe'])


def test_exceptions():

    with pytest.raises(AttributeError, match=r'Please use the FileSource class'):
        s = openmc.IndependentSource()
        s.file = 'my_file'

    with pytest.raises(AttributeError, match=r'Please use the CompiledSource class'):
        s = openmc.IndependentSource()
        s.library = 'my_library'

    with pytest.raises(AttributeError, match=r'Please use the CompiledSource class'):
        s = openmc.IndependentSource()
        s.parameters = 'my_params'

    with pytest.warns(FutureWarning, match=r'in favor of \'IndependentSource\''):
        s = openmc.Source()

    with pytest.raises(AttributeError, match=r'has no attribute \'frisbee\''):
        s = openmc.IndependentSource()
        s.frisbee
