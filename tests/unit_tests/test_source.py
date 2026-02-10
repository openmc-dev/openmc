from collections import Counter
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
    positions = [(1, 0, 2), (0, 1, 0), (0, 0, 3), (4, 9, 2)]
    strengths = [1, 2, 3, 4]

    space = openmc.stats.PointCloud(positions, strengths)
    np.testing.assert_equal(space.positions, positions)
    np.testing.assert_equal(space.strengths, strengths)

    src = openmc.IndependentSource(space=space)
    assert src.space == space
    np.testing.assert_equal(src.space.positions, positions)
    np.testing.assert_equal(src.space.strengths, strengths)

    elem = src.to_xml_element()
    src = openmc.IndependentSource.from_xml_element(elem)
    np.testing.assert_equal(src.space.positions, positions)
    np.testing.assert_equal(src.space.strengths, strengths)


def test_point_cloud_invalid():
    with pytest.raises(ValueError, match='2D'):
        openmc.stats.PointCloud([1, 0, 2, 0, 1, 0])

    with pytest.raises(ValueError, match='3 values'):
        openmc.stats.PointCloud([(1, 0, 2, 3), (4, 5, 2, 3)])

    with pytest.raises(ValueError, match='1D'):
        openmc.stats.PointCloud([(1, 0, 2), (4, 5, 2)], [(1, 2), (3, 4)])

    with pytest.raises(ValueError, match='same length'):
        openmc.stats.PointCloud([(1, 0, 2), (4, 5, 2)], [1, 2, 4])


def test_point_cloud_strengths(run_in_tmpdir, sphere_box_model):
    positions = [(1., 0., 2.), (0., 1., 0.), (0., 0., 3.), (-1., -1., 2.)]
    strengths = [1, 2, 3, 4]
    space = openmc.stats.PointCloud(positions, strengths)

    model = sphere_box_model[0]
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=space)

    try:
        model.init_lib()
        n_samples = 50_000
        sites = openmc.lib.sample_external_source(n_samples)
    finally:
        model.finalize_lib()

    count = Counter(s.r for s in sites)
    for i, (strength, position) in enumerate(zip(strengths, positions)):
        sampled_strength = count[position] / n_samples
        expected_strength = pytest.approx(strength/sum(strengths), abs=0.02)
        assert sampled_strength == expected_strength, f'Strength incorrect for {positions[i]}'


def test_source_file():
    filename = 'source.h5'
    src = openmc.FileSource(path=filename)
    assert src.path.name == filename

    elem = src.to_xml_element()
    assert 'strength' in elem.attrib
    assert 'file' in elem.attrib


def test_source_dlopen():
    library = 'libsource.so'
    src = openmc.CompiledSource(library)
    assert src.library.name == library

    elem = src.to_xml_element()
    assert 'library' in elem.attrib


def test_coincident_source_xml_roundtrip():
    # Create a coincident source (e.g., Co-60 two-gamma emission)
    gamma1 = openmc.IndependentSource(
        energy=openmc.stats.Discrete([1.17e6], [1.0]),
        angle=openmc.stats.Isotropic(),
        particle='photon'
    )
    gamma2 = openmc.IndependentSource(
        energy=openmc.stats.Discrete([1.33e6], [1.0]),
        angle=openmc.stats.Isotropic(),
        particle='photon'
    )
    space = openmc.stats.Point((1.0, 2.0, 3.0))
    time = openmc.stats.Discrete([0.0], [1.0])
    src = openmc.CoincidentSource(
        space=space,
        time=time,
        sources=[gamma1, gamma2],
        strength=5.0
    )

    # Write to XML element
    elem = src.to_xml_element()

    # Read back from XML element
    new_src = openmc.SourceBase.from_xml_element(elem)
    assert isinstance(new_src, openmc.CoincidentSource)
    assert new_src.type == 'coincident'
    assert new_src.strength == approx(5.0)

    # Check spatial distribution preserved
    assert isinstance(new_src.space, openmc.stats.Point)
    np.testing.assert_allclose(new_src.space.xyz, [1.0, 2.0, 3.0])

    # Check time distribution preserved
    assert isinstance(new_src.time, openmc.stats.Discrete)
    np.testing.assert_allclose(new_src.time.x, [0.0])
    np.testing.assert_allclose(new_src.time.p, [1.0])

    # Check sub-sources preserved
    assert len(new_src.sources) == 2
    for i, (orig, new) in enumerate(zip([gamma1, gamma2], new_src.sources)):
        assert isinstance(new, openmc.IndependentSource)
        assert new.particle == orig.particle
        assert isinstance(new.energy, openmc.stats.Discrete)
        np.testing.assert_allclose(new.energy.x, orig.energy.x)
        np.testing.assert_allclose(new.energy.p, orig.energy.p)


def test_coincident_source_probabilities_roundtrip():
    # Create a coincident source with emission probabilities
    gamma1 = openmc.IndependentSource(
        energy=openmc.stats.Discrete([1.17e6], [1.0]),
        angle=openmc.stats.Isotropic(),
        particle='photon'
    )
    gamma2 = openmc.IndependentSource(
        energy=openmc.stats.Discrete([1.33e6], [1.0]),
        angle=openmc.stats.Isotropic(),
        particle='photon'
    )
    src = openmc.CoincidentSource(
        space=openmc.stats.Point((0, 0, 0)),
        sources=[gamma1, gamma2],
        probabilities=[1.0, 0.5]
    )

    # Write to XML and read back
    elem = src.to_xml_element()
    new_src = openmc.SourceBase.from_xml_element(elem)
    assert isinstance(new_src, openmc.CoincidentSource)
    assert new_src.probabilities == [1.0, 0.5]
    assert len(new_src.sources) == 2

    # Verify the probability attribute is only written when != 1.0
    sub_elems = list(elem.iterchildren('source'))
    assert sub_elems[0].get('probability') is None  # 1.0 omitted
    assert sub_elems[1].get('probability') == '0.5'


def test_coincident_source_probabilities_default():
    # Without probabilities, they should be None
    gamma1 = openmc.IndependentSource(particle='photon')
    gamma2 = openmc.IndependentSource(particle='photon')
    src = openmc.CoincidentSource(
        space=openmc.stats.Point(),
        sources=[gamma1, gamma2]
    )
    assert src.probabilities is None

    # XML roundtrip should also give None
    elem = src.to_xml_element()
    new_src = openmc.SourceBase.from_xml_element(elem)
    assert new_src.probabilities is None


def test_coincident_source_validation():
    # Must have at least 2 sub-sources
    gamma1 = openmc.IndependentSource(particle='photon')
    with pytest.raises(ValueError, match='at least 2'):
        openmc.CoincidentSource(sources=[gamma1])


def test_coincident_source_probabilities_validation():
    gamma1 = openmc.IndependentSource(particle='photon')
    gamma2 = openmc.IndependentSource(particle='photon')
    src = openmc.CoincidentSource(
        space=openmc.stats.Point(),
        sources=[gamma1, gamma2]
    )

    # Wrong length
    with pytest.raises(ValueError, match='Length of probabilities'):
        src.probabilities = [0.5]

    # Out of range (0)
    with pytest.raises(ValueError, match='not in the range'):
        src.probabilities = [0.0, 0.5]

    # Out of range (> 1)
    with pytest.raises(ValueError, match='not in the range'):
        src.probabilities = [1.5, 0.5]


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
def test_rejection_fraction(run_in_tmpdir):
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    w = 0.25
    rpp1 = openmc.model.RectangularParallelepiped(
        -w/2, w/2, -w/2, w/2, -w/2, w/2)
    rpp2 = openmc.model.RectangularParallelepiped(
        -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=mat, region=-rpp1)
    cell2 = openmc.Cell(region=+rpp1 & -rpp2)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell1, cell2])

    # Create a box source over a 1 cm³ volume that is constrained to the source
    # cell of volume (0.25 cm)³ = 0.0125 cm³, which means the default rejection
    # fraction of 0.05 won't work
    model.settings.particles = 1000
    model.settings.batches = 1
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(*(-rpp2).bounding_box),
        constraints={'domains': [cell1]}
    )
    with pytest.raises(RuntimeError, match='Too few source sites'):
        model.run(openmc_exec=config['exe'])

    # With a source rejection fraction below 0.0125, the simulation should run
    model.settings.source_rejection_fraction = 0.005
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
