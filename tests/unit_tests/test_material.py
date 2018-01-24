import openmc
import openmc.model
import openmc.stats
import pytest


@pytest.fixture(scope='module')
def uo2():
    m = openmc.Material(material_id=100, name='UO2')
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('O16', 2.0)
    m.set_density('g/cm3', 10.0)
    m.depletable = True
    return m


@pytest.fixture(scope='module')
def sphere_model():
    model = openmc.model.Model()

    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.set_density('g/cm3', 1.0)
    model.materials.append(m)

    sph = openmc.Sphere(boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-sph)
    model.geometry.root_universe = openmc.Universe(cells=[c])

    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.Source(space=openmc.stats.Point())
    ll, ur = c.region.bounding_box
    model.settings.volume_calculations = [
        openmc.VolumeCalculation(domains=[m], samples=1000,
                                 lower_left=ll, upper_right=ur)
    ]
    return model


@pytest.fixture
def run_in_tmpdir(tmpdir):
    orig = tmpdir.chdir()
    try:
        yield
    finally:
        orig.chdir()


def test_attributes(uo2):
    assert uo2.name == 'UO2'
    assert uo2.id == 100
    assert uo2.depletable


def test_nuclides(uo2):
    """Test adding/removing nuclides."""
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    with pytest.raises(ValueError):
        m.add_nuclide('H1', '1.0')
    with pytest.raises(ValueError):
        m.add_nuclide(1.0, 'H1')
    with pytest.raises(ValueError):
        m.add_nuclide('H1', 1.0, 'oa')
    m.remove_nuclide('U235')


def test_elements():
    """Test adding elements."""
    m = openmc.Material()
    m.add_element('Zr', 1.0)
    m.add_element('U', 1.0, enrichment=4.5)
    with pytest.raises(ValueError):
        m.add_element('U', 1.0, enrichment=100.0)
    with pytest.raises(ValueError):
        m.add_element('Pu', 1.0, enrichment=3.0)


def test_density():
    m = openmc.Material()
    for unit in ['g/cm3', 'g/cc', 'kg/m3', 'atom/b-cm', 'atom/cm3']:
        m.set_density(unit, 1.0)
    with pytest.raises(ValueError):
        m.set_density('g/litre', 1.0)


def test_salphabeta():
    m = openmc.Material()
    m.add_s_alpha_beta('c_H_in_H2O', 0.5)


def test_repr():
    m = openmc.Material()
    m.add_nuclide('Zr90', 1.0)
    m.add_nuclide('H2', 0.5)
    m.add_s_alpha_beta('c_D_in_D2O')
    m.set_density('sum')
    m.temperature = 600.0
    repr(m)


def test_macroscopic():
    m = openmc.Material(name='UO2')
    m.add_macroscopic('UO2')
    with pytest.raises(ValueError):
        m.add_nuclide('H1', 1.0)
    with pytest.raises(ValueError):
        m.add_element('O', 1.0)
    with pytest.raises(ValueError):
        m.add_macroscopic('Other')

    m2 = openmc.Material()
    m2.add_nuclide('He4', 1.0)
    with pytest.raises(ValueError):
        m2.add_macroscopic('UO2')

    m.remove_macroscopic('UO2')


def test_isotropic():
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('O16', 2.0)
    m.make_isotropic_in_lab()
    assert m.isotropic == ['U235', 'O16']
    m.isotropic = ['O16']
    assert m.isotropic == ['O16']


def test_volume(run_in_tmpdir, sphere_model):
    """Test adding volume information from a volume calculation."""
    sphere_model.export_to_xml()
    openmc.calculate_volumes()
    volume_calc = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
    sphere_model.materials[0].add_volume_information(volume_calc)


def test_get_nuclide_densities(uo2):
    nucs = uo2.get_nuclide_densities()
    for nuc, density, density_type in nucs.values():
        assert nuc in ('U235', 'O16')
        assert density > 0
        assert density_type in ('ao', 'wo')


def test_get_nuclide_atom_densities(uo2):
    nucs = uo2.get_nuclide_atom_densities()
    for nuc, density in nucs.values():
        assert nuc in ('U235', 'O16')
        assert density > 0


def test_materials(run_in_tmpdir):
    m1 = openmc.Material()
    m1.add_nuclide('U235', 1.0, 'wo')
    m1.add_nuclide('O16', 2.0, 'wo')
    m1.set_density('g/cm3', 10.0)
    m1.depletable = True
    m1.temperature = 900.0

    m2 = openmc.Material()
    m2.add_nuclide('H1', 2.0)
    m2.add_nuclide('O16', 1.0)
    m2.add_s_alpha_beta('c_H_in_H2O')
    m2.set_density('kg/m3', 1000.0)

    mats = openmc.Materials([m1, m2])
    mats.cross_sections = '/some/fake/cross_sections.xml'
    mats.multipole_library = '/some/awesome/mp_lib/'
    mats.export_to_xml()
