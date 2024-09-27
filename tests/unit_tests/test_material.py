from collections import defaultdict
from pathlib import Path

import pytest

import openmc
from openmc.data import decay_photon_energy
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

def test_add_components():
    """Test adding multipe elements or nuclides at once"""
    m = openmc.Material()
    components = {'H1': 2.0,
                  'O16': 1.0,
                  'Zr': 1.0,
                  'O': 1.0,
                  'Ag110_m1': 1.0,
                  'U': {'percent': 1.0,
                        'enrichment': 4.5},
                  'Li': {'percent': 1.0,
                         'enrichment': 60.0,
                         'enrichment_target': 'Li7'},
                  'H': {'percent': 1.0,
                        'enrichment': 50.0,
                        'enrichment_target': 'H2',
                        'enrichment_type': 'wo'}}
    m.add_components(components)
    with pytest.raises(ValueError):
        m.add_components({'U': {'percent': 1.0,
                                'enrichment': 100.0}})
    with pytest.raises(ValueError):
        m.add_components({'Pu': {'percent': 1.0,
                                 'enrichment': 3.0}})
    with pytest.raises(ValueError):
        m.add_components({'U': {'percent': 1.0,
                                'enrichment': 70.0,
                                'enrichment_target':'U235'}})
    with pytest.raises(ValueError):
        m.add_components({'He': {'percent': 1.0,
                                 'enrichment': 17.0,
                                 'enrichment_target': 'He6'}})
    with pytest.raises(ValueError):
        m.add_components({'li': 1.0})  # should fail as 1st char is lowercase
    with pytest.raises(ValueError):
        m.add_components({'LI': 1.0})  # should fail as 2nd char is uppercase
    with pytest.raises(ValueError):
        m.add_components({'Xx': 1.0})  # should fail as Xx is not an element
    with pytest.raises(ValueError):
        m.add_components({'n': 1.0})  # check to avoid n for neutron being accepted
    with pytest.raises(TypeError):
        m.add_components({'H1': '1.0'})
    with pytest.raises(TypeError):
        m.add_components({1.0: 'H1'}, percent_type = 'wo')
    with pytest.raises(ValueError):
        m.add_components({'H1': 1.0}, percent_type = 'oa')

def test_nuclides_to_ignore(run_in_tmpdir):
    """Test nuclides_to_ignore when exporting a material to XML"""
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    m.add_nuclide('H1', 1.0)
    m.add_nuclide('O16', 1.0)

    mats = openmc.Materials([m])
    mats.export_to_xml(nuclides_to_ignore=['H1'])

    test_mats = openmc.Materials.from_xml()
    assert 'H1' not in test_mats[0].get_nuclides()

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


def test_remove_elements():
    """Test removing elements."""
    m = openmc.Material()
    for elem, percent in [('Li', 1.0), ('Be', 1.0)]:
        m.add_element(elem, percent)
    m.remove_element('Li')
    assert len(m.nuclides) == 1
    assert m.nuclides[0].name == 'Be9'
    assert m.nuclides[0].percent == 1.0


def test_add_element():
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
    with pytest.raises(ValueError):
        m.add_element('li', 1.0)  # should fail as 1st char is lowercase
    with pytest.raises(ValueError):
        m.add_element('LI', 1.0)  # should fail as 2nd char is uppercase
    with pytest.raises(ValueError):
        m.add_element('Xx', 1.0)  # should fail as Xx is not an element
    with pytest.raises(ValueError):
        m.add_element('n', 1.0)  # check to avoid n for neutron being accepted

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
    for nuclide, adens in m.get_nuclide_atom_densities().items():
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
        assert nuc_dens[nuclide] == pytest.approx(ref_dens[nuclide], 1e-2)

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
        assert nuc_dens[nuclide] == pytest.approx(ref_dens[nuclide], 1e-2)

    # testing the use of brackets
    m = openmc.Material()
    m.add_elements_from_formula('Mg2(NO3)2')

    # checking the ratio of elements is 2:2:6 for Mg:N:O
    elem = defaultdict(float)
    for nuclide, adens in m.get_nuclide_atom_densities().items():
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
        assert nuc_dens[nuclide] == pytest.approx(ref_dens[nuclide], 1e-2)

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


def test_get_nuclides():
    mat = openmc.Material()

    mat.add_nuclide('Li6', 1.0)
    assert mat.get_nuclides() == ['Li6']
    assert mat.get_nuclides(element='Li') == ['Li6']
    assert mat.get_nuclides(element='Be') == []

    mat.add_element('Li', 1.0)
    assert mat.get_nuclides() == ['Li6', 'Li7']
    assert mat.get_nuclides(element='Be') == []

    mat.add_element('Be', 1.0)
    assert mat.get_nuclides() == ['Li6', 'Li7', 'Be9']
    assert mat.get_nuclides(element='Be') == ['Be9']


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
    for nuc, density in uo2.get_nuclide_atom_densities().items():
        assert nuc in ('U235', 'O16')
        assert density > 0


def test_get_nuclide_atom_densities_specific(uo2):
    one_nuc = uo2.get_nuclide_atom_densities(nuclide='O16')
    assert list(one_nuc.keys()) == ['O16']
    assert list(one_nuc.values())[0] > 0

    all_nuc = uo2.get_nuclide_atom_densities()
    assert all_nuc['O16'] == one_nuc['O16']


def test_get_nuclide_atoms():
    mat = openmc.Material()
    mat.add_nuclide('Li6', 1.0)
    mat.set_density('atom/cm3', 3.26e20)
    mat.volume = 100.0

    atoms = mat.get_nuclide_atoms()
    assert atoms['Li6'] == pytest.approx(mat.density * mat.volume)

    atoms = mat.get_nuclide_atoms(volume=10.0)
    assert atoms['Li6'] == pytest.approx(mat.density * 10.0)


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

    # Test with volume specified as argument
    assert m.get_mass('Zr90', volume=1.0) == pytest.approx(1.0)


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
        assert nuc_dens[nuclide] == pytest.approx(ref_dens[nuclide], 1e-2)
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


def test_get_activity():
    """Tests the activity of stable, metastable and active materials"""

    # Creates a material with stable isotopes to check the activity is 0
    m1 = openmc.Material()
    m1.add_element("Fe", 0.7)
    m1.add_element("Li", 0.3)
    m1.set_density('g/cm3', 1.5)
    # activity in Bq/cc and Bq/g should not require volume setting
    assert m1.get_activity(units='Bq/cm3') == 0
    assert m1.get_activity(units='Bq/g') == 0
    m1.volume = 1
    assert m1.get_activity(units='Bq') == 0

    # Checks that 1g of tritium has the correct activity scaling
    m2 = openmc.Material()
    m2.add_nuclide("H3", 1)
    m2.set_density('g/cm3', 1)
    m2.volume = 1
    assert pytest.approx(m2.get_activity(units='Bq')) == 3.559778e14
    m2.set_density('g/cm3', 2)
    assert pytest.approx(m2.get_activity(units='Bq')) == 3.559778e14*2
    m2.volume = 3
    assert pytest.approx(m2.get_activity(units='Bq')) == 3.559778e14*2*3

    # Checks that 1 mol of a metastable nuclides has the correct activity
    m3 = openmc.Material()
    m3.add_nuclide("Tc99_m1", 1)
    m3.set_density('g/cm3', 1)
    m3.volume = 98.9
    assert pytest.approx(m3.get_activity(units='Bq'), rel=0.001) == 1.93e19

    # Checks that specific and volumetric activity of tritium are correct
    m4 = openmc.Material()
    m4.add_nuclide("H3", 1)
    m4.set_density('g/cm3', 1.5)
    assert pytest.approx(m4.get_activity(units='Bq/g')) == 355978108155965.94  # [Bq/g]
    assert pytest.approx(m4.get_activity(units='Bq/g', by_nuclide=True)["H3"]) == 355978108155965.94  # [Bq/g]
    assert pytest.approx(m4.get_activity(units='Bq/cm3')) == 355978108155965.94*3/2 # [Bq/cc]
    assert pytest.approx(m4.get_activity(units='Bq/cm3', by_nuclide=True)["H3"]) == 355978108155965.94*3/2 # [Bq/cc]
    # volume is required to calculate total activity
    m4.volume = 10.
    assert pytest.approx(m4.get_activity(units='Bq')) == 355978108155965.94*3/2*10 # [Bq]

    # Test with volume specified as argument
    assert pytest.approx(m4.get_activity(units='Bq', volume=1.0)) == 355978108155965.94*3/2


def test_get_decay_heat():
    # Set chain file for testing
    openmc.config['chain_file'] = Path(__file__).parents[1] / 'chain_simple.xml'

    """Tests the decay heat of stable, metastable and active materials"""
    m1 = openmc.Material()
    m1.add_nuclide("U235", 0.2)
    m1.add_nuclide("U238", 0.8)
    m1.set_density('g/cm3', 10.5)
    # decay heat in W/cc and W/g should not require volume setting
    assert m1.get_decay_heat(units='W/cm3') == 0
    assert m1.get_decay_heat(units='W/g') == 0
    m1.volume = 1
    assert m1.get_decay_heat(units='W') == 0

    # Checks that 1g of tritium has the correct decay heat scaling
    m2 = openmc.Material()
    m2.add_nuclide("I135", 1)
    m2.set_density('g/cm3', 1)
    m2.volume = 1
    assert pytest.approx(m2.get_decay_heat(units='W')) == 40175.15720273193
    m2.set_density('g/cm3', 2)
    assert pytest.approx(m2.get_decay_heat(units='W')) == 40175.15720273193*2
    m2.volume = 3
    assert pytest.approx(m2.get_decay_heat(units='W')) == 40175.15720273193*2*3

    # Checks that 1 mol of a metastable nuclides has the correct decay heat
    m3 = openmc.Material()
    m3.add_nuclide("Xe135", 1)
    m3.set_density('g/cm3', 1)
    m3.volume = 98.9
    assert pytest.approx(m3.get_decay_heat(units='W'), rel=0.001) == 846181.2921143445

    # Checks that specific and volumetric decay heat of tritium are correct
    m4 = openmc.Material()
    m4.add_nuclide("I135", 1)
    m4.set_density('g/cm3', 1.5)
    assert pytest.approx(m4.get_decay_heat(units='W/g')) == 40175.15720273193 # [W/g]
    assert pytest.approx(m4.get_decay_heat(units='W/g', by_nuclide=True)["I135"]) == 40175.15720273193 # [W/g]
    assert pytest.approx(m4.get_decay_heat(units='W/cm3')) == 40175.15720273193*3/2 # [W/cc]
    assert pytest.approx(m4.get_decay_heat(units='W/cm3', by_nuclide=True)["I135"]) == 40175.15720273193*3/2 #[W/cc]
    # volume is required to calculate total decay heat
    m4.volume = 10.
    assert pytest.approx(m4.get_decay_heat(units='W')) == 40175.15720273193*3/2*10 # [W]

    # Test with volume specified as argument
    assert pytest.approx(m4.get_decay_heat(units='W', volume=1.0)) == 40175.15720273193*3/2


def test_decay_photon_energy():
    # Set chain file for testing
    openmc.config['chain_file'] = Path(__file__).parents[1] / 'chain_simple.xml'

    # Material representing single atom of I135 and Cs135
    m = openmc.Material()
    m.add_nuclide('I135', 1.0e-24)
    m.add_nuclide('Cs135', 1.0e-24)
    m.volume = 1.0

    # Get decay photon source and make sure it's the right type
    src = m.get_decay_photon_energy()
    assert isinstance(src, openmc.stats.Discrete)

    # Make sure units/volume work as expected
    src_v2 = m.get_decay_photon_energy(volume=2.0)
    assert src.p * 2.0 == pytest.approx(src_v2.p)
    src_per_cm3 = m.get_decay_photon_energy(units='Bq/cm3', volume=100.0)
    assert (src.p == src_per_cm3.p).all()

    # If we add Xe135 (which has a tabular distribution), the photon source
    # should be a mixture distribution
    m.add_nuclide('Xe135', 1.0e-24)
    src = m.get_decay_photon_energy()
    assert isinstance(src, openmc.stats.Mixture)

    # With a single atom of each, the intensity of the photon source should be
    # equal to the sum of the intensities for each nuclide
    def intensity(src):
        return src.integral() if src is not None else 0.0

    assert src.integral() == pytest.approx(sum(
        intensity(decay_photon_energy(nuc)) for nuc in m.get_nuclides()
    ), rel=1e-3)

    # When the clipping threshold is zero, the intensities should match exactly
    src = m.get_decay_photon_energy(0.0)
    assert src.integral() == pytest.approx(sum(
        intensity(decay_photon_energy(nuc)) for nuc in m.get_nuclides()
    ))

    # A material with no unstable nuclides should have no decay photon source
    stable = openmc.Material()
    stable.add_nuclide('Gd156', 1.0)
    stable.volume = 1.0
    assert stable.get_decay_photon_energy() is None


def test_avoid_subnormal(run_in_tmpdir):
    # Write a materials.xml with a material that has a nuclide density that is
    # represented as a subnormal floating point value
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.add_nuclide('H2', 1.0e-315)
    mats = openmc.Materials([mat])
    mats.export_to_xml()

    # When read back in, the density should be zero
    mats = openmc.Materials.from_xml()
    assert mats[0].get_nuclide_atom_densities()['H2'] == 0.0
