from importlib import resources
import json

import lxml.etree as ET
import pytest

import openmc
import openmc.material as material_module


@pytest.fixture
def material_library_registry():
    """Restore the material library registry and cache after a test."""
    registry = material_module._MATERIAL_LIBRARIES.copy()
    try:
        yield
    finally:
        material_module._MATERIAL_LIBRARIES.clear()
        material_module._MATERIAL_LIBRARIES.update(registry)
        material_module._load_material_library.cache_clear()


@pytest.fixture
def cross_sections(tmp_path):
    """Cross section index covering the natural elements used in tests."""
    root = ET.Element('cross_sections')
    for nuclide in (
        'H1', 'H2', 'O16', 'O17', 'Na23', 'Si28', 'Si29', 'Si30', 'Cl35',
        'Cl37', 'Y89', 'Ce136', 'Ce138', 'Ce140', 'Ce142', 'Cs133',
        'Lu175', 'Lu176'
    ):
        ET.SubElement(
            root,
            'library',
            materials=nuclide,
            path=f'{nuclide}.h5',
            type='neutron',
        )
    path = tmp_path / 'cross_sections.xml'
    ET.ElementTree(root).write(path)
    return path


def test_natural_elements(cross_sections):
    """Natural elements use nuclides present in the cross section library."""
    with openmc.config.patch('cross_sections', cross_sections):
        material = openmc.Material.from_library('Sodium Oxide')

    assert material.name == 'Sodium Oxide'
    assert material.density == pytest.approx(2.27)
    assert material.density_units == 'g/cm3'

    nuclides = {nuc.name: nuc.percent for nuc in material.nuclides}
    assert set(nuclides) == {'O16', 'O17', 'Na23'}
    assert nuclides['O16'] + nuclides['O17'] == pytest.approx(0.333333)
    assert nuclides['Na23'] == pytest.approx(0.666667)


def test_isotopic_components(cross_sections):
    """Non-natural isotopic compositions are retained explicitly."""
    he3 = openmc.Material.from_library('He-3 Proportional Gas')
    assert he3.nuclides == [('He3', 1.0, 'ao')]

    with openmc.config.patch('cross_sections', cross_sections):
        clyc = openmc.Material.from_library(
            'Cesium Lithium Yttrium Chloride (CLYC) with 95% Li6 Enrichment'
        )
    lithium = {
        nuc.name: nuc.percent
        for nuc in clyc.nuclides
        if nuc.name.startswith('Li')
    }
    assert lithium == pytest.approx({'Li6': 0.095682, 'Li7': 0.004318})

    leu = openmc.Material.from_library('Uranium, Low Enriched (LEU)')
    uranium = {nuc.name: nuc.percent for nuc in leu.nuclides}
    assert uranium == pytest.approx({
        'U234': 0.000271,
        'U235': 0.030372,
        'U236': 0.000139,
        'U238': 0.969217,
    })

    aged_pu = openmc.Material.from_library(
        'Plutonium, Aged WGPu (A: 4-7% Pu240)'
    )
    assert 'Am241' in aged_pu.get_nuclides()
    assert sum(nuc.percent for nuc in aged_pu.nuclides) == pytest.approx(1.0)


def test_independent_materials():
    """Each lookup creates an independent material with a fresh ID."""
    first = openmc.Material.from_library('He-3 Proportional Gas')
    second = openmc.Material.from_library('He-3 Proportional Gas')

    assert first is not second
    assert first.id != second.id
    first.add_nuclide('H1', 1.0)
    assert second.get_nuclides() == ['He3']


def test_constructor_kwargs():
    """Keyword arguments are forwarded to the Material constructor."""
    material = openmc.Material.from_library(
        'He-3 Proportional Gas',
        material_id=987654,
        name='Helium-3 detector gas',
        temperature=293.15,
        volume=10.0,
        depletable=True,
    )

    assert material.id == 987654
    assert material.name == 'Helium-3 detector gas'
    assert material.temperature == 293.15
    assert material.volume == 10.0
    assert material.depletable
    assert material.density == pytest.approx(0.0001252645124733361)
    assert material.nuclides == [('He3', 1.0, 'ao')]


def test_library_material_names():
    """Library names are sorted and do not consume a material ID."""
    next_id = openmc.Material.next_id
    names = openmc.Material.get_library_material_names()

    assert isinstance(names, tuple)
    assert names == tuple(sorted(names))
    assert len(names) == 411
    assert 'Lutetium Yttrium OxyorthoSilicate: 0.5 atom% Cerium (LYSO)' in names
    assert 'Glass Scintillator, Li Doped  (GS1, GS2, GS3)' not in names
    assert 'Glass Scintillator, Li Doped (GS1, GS2, GS3)' in names
    assert 'Sodium Oxide' in names
    assert openmc.Material.next_id == next_id


def test_library_name_whitespace():
    """Whitespace differences do not prevent material lookup."""
    stored_name = 'Glass Scintillator, Li Doped (GS1, GS2, GS3)'
    query_name = 'Glass Scintillator, Li Doped  (GS1, GS2, GS3)'
    stored = openmc.Material.from_library(stored_name)
    queried = openmc.Material.from_library(query_name)

    assert queried.name == query_name
    assert queried.density == stored.density
    assert queried.get_nuclides() == stored.get_nuclides()


def test_register_library(
    tmp_path, cross_sections, material_library_registry
):
    """A custom library can be registered, listed, and loaded."""
    path = tmp_path / 'custom_materials.json'
    path.write_text(json.dumps({
        'schema_version': 1,
        'density_units': 'g/cm3',
        'percent_type': 'ao',
        'materials': {
            'Custom  Water': {
                'density': 0.95,
                'elements': {'H': 0.666667, 'O': 0.333333},
            },
        },
    }))

    openmc.Material.register_library('custom', path)
    assert openmc.Material.get_library_material_names('custom') == (
        'Custom Water',
    )

    with openmc.config.patch('cross_sections', cross_sections):
        material = openmc.Material.from_library(
            'Custom Water',
            library='custom',
            material_id=987653,
            temperature=600.0,
        )

    assert material.id == 987653
    assert material.temperature == 600.0
    assert material.density == pytest.approx(0.95)
    fractions = {nuc.name: nuc.percent for nuc in material.nuclides}
    assert fractions['H1'] + fractions['H2'] == pytest.approx(0.666667)
    assert fractions['O16'] + fractions['O17'] == pytest.approx(0.333333)

    with pytest.raises(ValueError, match='already registered'):
        openmc.Material.register_library('custom', path)
    assert openmc.Material.get_library_material_names('custom') == (
        'Custom Water',
    )


def test_register_invalid_library(tmp_path, material_library_registry):
    """Registration validates names, files, and top-level schemas."""
    valid_data = {
        'schema_version': 1,
        'density_units': 'g/cm3',
        'percent_type': 'ao',
        'materials': {
            'Hydrogen': {
                'density': 0.1,
                'nuclides': {'H1': 1.0},
            },
        },
    }

    valid_path = tmp_path / 'valid.json'
    valid_path.write_text(json.dumps(valid_data))
    with pytest.raises(ValueError, match='cannot be empty'):
        openmc.Material.register_library('', valid_path)
    with pytest.raises(ValueError, match='already registered'):
        openmc.Material.register_library('pnnl_v2', valid_path)

    with pytest.raises(RuntimeError, match='Could not load'):
        openmc.Material.register_library('missing', tmp_path / 'missing.json')

    malformed_path = tmp_path / 'malformed.json'
    malformed_path.write_text('{')
    with pytest.raises(RuntimeError, match='Could not load'):
        openmc.Material.register_library('malformed', malformed_path)

    schema_path = tmp_path / 'schema.json'
    schema_path.write_text(json.dumps({**valid_data, 'schema_version': 2}))
    with pytest.raises(RuntimeError, match='unsupported schema version'):
        openmc.Material.register_library('schema', schema_path)

    invalid_path = tmp_path / 'invalid.json'
    invalid_data = {
        **valid_data,
        'materials': {
            'Invalid element': {
                'density': 0.1,
                'elements': {'Xx': 1.0},
            },
        },
    }
    invalid_path.write_text(json.dumps(invalid_data))
    openmc.Material.register_library('invalid', invalid_path)
    with pytest.raises(ValueError, match='not recognised'):
        openmc.Material.from_library('Invalid element', library='invalid')

    with pytest.raises(ValueError, match="Unknown material library 'missing'"):
        openmc.Material.get_library_material_names('missing')


def test_unknown_library_and_material():
    with pytest.raises(ValueError, match="Unknown material library 'unknown'"):
        openmc.Material.from_library('Sodium Oxide', library='unknown')

    with pytest.raises(ValueError, match="Material 'unknown' not found"):
        openmc.Material.from_library('unknown')


def test_pnnl_data():
    """Check invariants of the complete bundled PNNL material library."""
    path = resources.files('openmc.data').joinpath(
        'material_libraries/pnnl_v2.json'
    )
    data = json.loads(path.read_text(encoding='utf-8'))

    assert data['schema_version'] == 1
    assert data['density_units'] == 'g/cm3'
    assert data['percent_type'] == 'ao'
    assert data['source']['data_origin'] == (
        'PNNL Materials Compendium downloadable JSON'
    )
    assert 'LYSO' in data['source']['note']
    assert data['source']['report_sha256'] == (
        '72b26dba2c3b5583b86fe5d5fe27a43d2890331d0515ce787f1c18fd7321cee6'
    )
    assert data['source']['data_sha256'] == (
        '5db6f9ca58793659e73c66cfd624736ffa59846239cd423c11eb6dbe5a56607b'
    )
    lyso = data['materials'][
        'Lutetium Yttrium OxyorthoSilicate: 0.5 atom% Cerium (LYSO)'
    ]
    assert lyso['density'] == 7.25
    assert lyso['elements'] == {
        'O': 0.621875,
        'Si': 0.124375,
        'Y': 0.012438,
        'Lu': 0.236313,
        'Ce': 0.005,
    }
    assert 'PNNL-15870 Rev. 2' in lyso['source_note']
    assert len(data['materials']) == 411
    assert sum('nuclides' in mat for mat in data['materials'].values()) == 45

    element_symbols = set(openmc.data.ATOMIC_SYMBOL.values())
    for name, material in data['materials'].items():
        assert name
        assert material['density'] > 0.0

        elements = material.get('elements', {})
        nuclides = material.get('nuclides', {})
        assert elements or nuclides
        assert set(elements) <= element_symbols
        for nuclide in nuclides:
            openmc.data.zam(nuclide)

        fractions = [*elements.values(), *nuclides.values()]
        assert all(fraction > 0.0 for fraction in fractions)
        assert sum(fractions) == pytest.approx(1.0, abs=5.0e-6)
