from collections import defaultdict

import pytest

import openmc
import openmc.examples
import openmc.model
import openmc.stats


def test_attributes(uo2):
    assert uo2.name == 'UO2'
    assert uo2.id == 100
    assert uo2.depletable


def test_add_nuclide():
    """Test adding nuclides."""
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    with pytest.raises(TypeError):
        m.add_nuclide('H1', '1.0')
    with pytest.raises(TypeError):
        m.add_nuclide(1.0, 'H1')
    with pytest.raises(ValueError):
        m.add_nuclide('H1', 1.0, 'oa')


def test_remove_nuclide():
    """Test removing nuclides."""
    m = openmc.Material()
    for nuc, percent in [('H1', 1.0), ('H2', 1.0), ('H1', 2.0), ('H2', 2.0)]:
        m.add_nuclide(nuc, percent)
    m.remove_nuclide('H1')
    assert len(m.nuclides) == 2
    assert all(nuc.name == 'H2' for nuc in m.nuclides)
    assert m.nuclides[0].percent == 1.0
    assert m.nuclides[1].percent == 2.0


def test_elements():
    """Test adding elements."""
    m = openmc.Material()
    m.add_element('Zr', 1.0)
    m.add_element('U', 1.0, enrichment=4.5)
    m.add_element('Li', 1.0, enrichment=60.0, enrichment_target='Li7')
    m.add_element('H', 1.0, enrichment=50.0, enrichment_target='H2',
                  enrichment_type='wo')
    with pytest.raises(ValueError):
        m.add_element('U', 1.0, enrichment=100.0)
    with pytest.raises(ValueError):
        m.add_element('Pu', 1.0, enrichment=3.0)
    with pytest.raises(ValueError):
        m.add_element('U', 1.0, enrichment=70.0, enrichment_target='U235')
    with pytest.raises(ValueError):
        m.add_element('He', 1.0, enrichment=17.0, enrichment_target='He6')

def test_elements_by_name():
    """Test adding elements by name"""
    m = openmc.Material()
    m.add_element('woLfrAm', 1.0)
    with pytest.raises(ValueError):
        m.add_element('uranum', 1.0)
    m.add_element('uRaNiUm', 1.0)
    m.add_element('Aluminium', 1.0)
    a = openmc.Material()
    b = openmc.Material()
    c = openmc.Material()
    a.add_element('sulfur', 1.0)
    b.add_element('SulPhUR', 1.0)
    c.add_element('S', 1.0)
    assert a._nuclides == b._nuclides
    assert b._nuclides == c._nuclides


def test_add_elements_by_formula():
    """Test adding elements from a formula"""
    # testing the correct nuclides and elements are added to a material
    m = openmc.Material()
    m.add_elements_from_formula('Li4SiO4')
    # checking the ratio of elements is 4:1:4 for Li:Si:O
    elem = defaultdict(float)
    for nuclide, adens in m.get_nuclide_atom_densities().values():
        if nuclide.startswith("Li"):
            elem["Li"] += adens
        if nuclide.startswith("Si"):
            elem["Si"] += adens
        if nuclide.startswith("O"):
            elem["O"] += adens
    total_number_of_atoms = 9
    assert elem["Li"] == pytest.approx(4./total_number_of_atoms)
    assert elem["Si"] == pytest.approx(1./total_number_of_atoms)
    assert elem["O"] == pytest.approx(4/total_number_of_atoms)
    # testing the correct nuclides are added to the Material
    ref_dens = {'Li6': 0.033728, 'Li7': 0.410715,
                'Si28': 0.102477, 'Si29': 0.0052035, 'Si30': 0.0034301,
                'O16': 0.443386, 'O17': 0.000168}
    nuc_dens = m.get_nuclide_atom_densities()
    for nuclide in ref_dens:
        assert nuc_dens[nuclide][1] == pytest.approx(ref_dens[nuclide], 1e-2)

    # testing the correct nuclides are added to the Material when enriched
    m = openmc.Material()
    m.add_elements_from_formula('Li4SiO4',
                                enrichment=60.,
                                enrichment_target='Li6')
    ref_dens = {'Li6': 0.2666, 'Li7': 0.1777,
                'Si28': 0.102477, 'Si29': 0.0052035, 'Si30': 0.0034301,
                'O16': 0.443386, 'O17': 0.000168}
    nuc_dens = m.get_nuclide_atom_densities()
    for nuclide in ref_dens:
        assert nuc_dens[nuclide][1] == pytest.approx(ref_dens[nuclide], 1e-2)

    # testing the use of brackets
    m = openmc.Material()
    m.add_elements_from_formula('Mg2(NO3)2')

    # checking the ratio of elements is 2:2:6 for Mg:N:O
    elem = defaultdict(float)
    for nuclide, adens in m.get_nuclide_atom_densities().values():
        if nuclide.startswith("Mg"):
            elem["Mg"] += adens
        if nuclide.startswith("N"):
            elem["N"] += adens
        if nuclide.startswith("O"):
            elem["O"] += adens
    total_number_of_atoms = 10
    assert elem["Mg"] == pytest.approx(2./total_number_of_atoms)
    assert elem["N"] == pytest.approx(2./total_number_of_atoms)
    assert elem["O"] == pytest.approx(6/total_number_of_atoms)

    # testing the correct nuclides are added when brackets are used
    ref_dens = {'Mg24': 0.157902, 'Mg25': 0.02004, 'Mg26': 0.022058,
                'N14': 0.199267, 'N15': 0.000732,
                'O16': 0.599772, 'O17': 0.000227}
    nuc_dens = m.get_nuclide_atom_densities()
    for nuclide in ref_dens:
        assert nuc_dens[nuclide][1] == pytest.approx(ref_dens[nuclide], 1e-2)

    # testing non integer multiplier results in a value error
    m = openmc.Material()
    with pytest.raises(ValueError):
        m.add_elements_from_formula('Li4.2SiO4')

    # testing lowercase elements results in a value error
    m = openmc.Material()
    with pytest.raises(ValueError):
        m.add_elements_from_formula('li4SiO4')

    # testing lowercase elements results in a value error
    m = openmc.Material()
    with pytest.raises(ValueError):
        m.add_elements_from_formula('Li4Sio4')

    # testing incorrect character in formula results in a value error
    m = openmc.Material()
    with pytest.raises(ValueError):
        m.add_elements_from_formula('Li4$SiO4')

    # testing unequal opening and closing brackets
    m = openmc.Material()
    with pytest.raises(ValueError):
        m.add_elements_from_formula('Fe(H2O)4(OH)2)')


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


def test_macroscopic(run_in_tmpdir):
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

    # Make sure we can remove/add macroscopic
    m.remove_macroscopic('UO2')
    m.add_macroscopic('UO2')
    repr(m)

    # Make sure we can export a material with macroscopic data
    mats = openmc.Materials([m])
    mats.export_to_xml()


def test_paths():
    model = openmc.examples.pwr_assembly()
    model.geometry.determine_paths()
    fuel = model.materials[0]
    assert fuel.num_instances == 264
    assert len(fuel.paths) == 264


def test_isotropic():
    m1 = openmc.Material()
    m1.add_nuclide('U235', 1.0)
    m1.add_nuclide('O16', 2.0)
    m1.isotropic = ['O16']
    assert m1.isotropic == ['O16']

    m2 = openmc.Material()
    m2.add_nuclide('H1', 1.0)
    mats = openmc.Materials([m1, m2])
    mats.make_isotropic_in_lab()
    assert m1.isotropic == ['U235', 'O16']
    assert m2.isotropic == ['H1']


def test_get_elements():
    # test that zero elements exist on creation
    m = openmc.Material()
    assert len(m.get_elements()) == 0

    # test addition of a single element
    m.add_element('Li', 0.2)
    assert m.get_elements() == ["Li"]

    # test that adding the same element
    m.add_element('Li', 0.3)
    assert m.get_elements() == ["Li"]

    # test adding another element
    m.add_element('Si', 0.3)
    assert m.get_elements() == ["Li", "Si"]

    # test adding a third element
    m.add_element('O', 0.4)
    assert m.get_elements() == ["Li", "O", "Si"]
    # test removal of nuclides
    m.remove_nuclide('O16')
    m.remove_nuclide('O17')
    assert m.get_elements() == ["Li", "Si"]


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


def test_mass():
    m = openmc.Material()
    m.add_nuclide('Zr90', 1.0, 'wo')
    m.add_nuclide('U235', 1.0, 'wo')
    m.set_density('g/cm3', 2.0)
    m.volume = 10.0

    assert m.get_mass_density('Zr90') == pytest.approx(1.0)
    assert m.get_mass_density('U235') == pytest.approx(1.0)
    assert m.get_mass_density() == pytest.approx(2.0)

    assert m.get_mass('Zr90') == pytest.approx(10.0)
    assert m.get_mass('U235') == pytest.approx(10.0)
    assert m.get_mass() == pytest.approx(20.0)
    assert m.fissionable_mass == pytest.approx(10.0)


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
    mats.export_to_xml()


def test_borated_water():
    # Test against reference values from the BEAVRS benchmark.
    m = openmc.model.borated_water(975, 566.5, 15.51, material_id=50)
    assert m.density == pytest.approx(0.7405, 1e-3)
    assert m.temperature == pytest.approx(566.5)
    assert m._sab[0][0] == 'c_H_in_H2O'
    ref_dens = {'B10':8.0023e-06, 'B11':3.2210e-05, 'H1':4.9458e-02,
                'O16':2.4672e-02}
    nuc_dens = m.get_nuclide_atom_densities()
    for nuclide in ref_dens:
        assert nuc_dens[nuclide][1] == pytest.approx(ref_dens[nuclide], 1e-2)
    assert m.id == 50

    # Test the Celsius conversion.
    m = openmc.model.borated_water(975, 293.35, 15.51, 'C')
    assert m.density == pytest.approx(0.7405, 1e-3)

    # Test Fahrenheit and psi conversions.
    m = openmc.model.borated_water(975, 560.0, 2250.0, 'F', 'psi')
    assert m.density == pytest.approx(0.7405, 1e-3)

    # Test the density override
    m = openmc.model.borated_water(975, 566.5, 15.51, density=0.9)
    assert m.density == pytest.approx(0.9, 1e-3)


def test_from_xml(run_in_tmpdir):
    # Create a materials.xml file
    m1 = openmc.Material(1, 'water')
    m1.add_nuclide('H1', 1.0)
    m1.add_nuclide('O16', 2.0)
    m1.add_s_alpha_beta('c_H_in_H2O')
    m1.temperature = 300
    m1.volume = 100
    m1.set_density('g/cm3', 0.9)
    m1.isotropic = ['H1']
    m2 = openmc.Material(2, 'zirc')
    m2.add_nuclide('Zr90', 1.0, 'wo')
    m2.set_density('kg/m3', 10.0)
    m3 = openmc.Material(3)
    m3.add_nuclide('N14', 0.02)

    mats = openmc.Materials([m1, m2, m3])
    mats.cross_sections = 'fake_path.xml'
    mats.export_to_xml()

    # Regenerate materials from XML
    mats = openmc.Materials.from_xml()
    assert len(mats) == 3
    m1 = mats[0]
    assert m1.id == 1
    assert m1.name == 'water'
    assert m1.nuclides == [('H1', 1.0, 'ao'), ('O16', 2.0, 'ao')]
    assert m1.isotropic == ['H1']
    assert m1.temperature == 300
    assert m1.volume == 100
    m2 = mats[1]
    assert m2.nuclides == [('Zr90', 1.0, 'wo')]
    assert m2.density == 10.0
    assert m2.density_units == 'kg/m3'
    assert mats[2].density_units == 'sum'


def test_mix_materials():
    m1 = openmc.Material()
    m1.add_nuclide('U235', 1.)
    m1dens = 10.0
    m1amm = m1.average_molar_mass
    m1.set_density('g/cm3', m1dens)
    m2 = openmc.Material()
    m2.add_nuclide('Zr90', 1.)
    m2dens = 2.0
    m2amm = m2.average_molar_mass
    m2.set_density('g/cm3', m2dens)
    f0, f1 = 0.6, 0.4
    dens3 = (f0*m1amm + f1*m2amm) / (f0*m1amm/m1dens + f1*m2amm/m2dens)
    dens4 = 1. / (f0 / m1dens + f1 / m2dens)
    dens5 = f0*m1dens + f1*m2dens
    m3 = openmc.Material.mix_materials([m1, m2], [f0, f1], percent_type='ao')
    m4 = openmc.Material.mix_materials([m1, m2], [f0, f1], percent_type='wo')
    m5 = openmc.Material.mix_materials([m1, m2], [f0, f1], percent_type='vo')
    assert m3.density == pytest.approx(dens3)
    assert m4.density == pytest.approx(dens4)
    assert m5.density == pytest.approx(dens5)
