""" Tests for ExternalSourceRates class """

from pathlib import Path
from math import exp

import pytest
import numpy as np
import re

import openmc
from openmc.deplete import CoupledOperator
from openmc.deplete.transfer_rates import ExternalSourceRates
from openmc.deplete.abc import (_SECONDS_PER_MINUTE, _SECONDS_PER_HOUR,
                                _SECONDS_PER_DAY, _SECONDS_PER_JULIAN_YEAR)

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    f = openmc.Material(name="f")
    f.add_element("U", 1, percent_type="ao", enrichment=4.25)
    f.add_element("O", 2)
    f.set_density("g/cc", 10.4)

    w = openmc.Material(name="w")
    w.add_element("O", 1)
    w.add_element("H", 2)
    w.set_density("g/cc", 1.0)
    w.depletable = True

    # material just to test multiple destination material
    h = openmc.Material(name="h")
    h.add_element("He", 1)
    h.set_density("g/cc", 1.78e-4)

    radii = [0.42, 0.45]
    f.volume = np.pi * radii[0] ** 2
    w.volume = np.pi * (radii[1]**2 - radii[0]**2)
    h.volume = 1
    materials = openmc.Materials([f, w, h])

    surf_f = openmc.Sphere(r=radii[0])
    surf_w = openmc.Sphere(r=radii[1], boundary_type='vacuum')
    surf_h = openmc.Sphere(x0=10, r=1, boundary_type='vacuum')
    cell_f = openmc.Cell(fill=f, region=-surf_f)
    cell_w = openmc.Cell(fill=w, region=+surf_f & -surf_w)
    cell_h = openmc.Cell(fill=h, region=-surf_h)
    geometry = openmc.Geometry([cell_f, cell_w, cell_h])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 10
    settings.batches = 50

    return openmc.Model(geometry, materials, settings)

@pytest.mark.parametrize(
"case_name, external_source_vectors, external_source_rate, timesteps", [
    ('elements', [{'U': 0.9, 'Xe': 0.1}], 1, None),
    ('nuclides', [{'I135': 0.1, 'Gd156': 0.3, 'Gd157': 0.6}], 1, None),
    ('nuclides_elements', [{'I135': 0.01, 'Gd156': 0.1, 'Gd157': 0.01, 'U': 0.8,
                           'Xe': 0.08}], 1, None),
    ('elements_nuclides', [{'U': 0.78, 'Xe': 0.1, 'I135': 0.01, 'Gd156': 0.1,
                           'Gd157': 0.01}], 1, None),
    ('multiple_vectors', [{'U': 1.}, {'Xe': 1}], 1, None),
    ('timesteps', [{'U': 0.9, 'Xe': 0.1}], 1, [1]),
    ('rates_invalid_1', [{'Gb': 1.}], 1, None),
    ('rates_invalid_2', [{'Pu': 1.}], 1, None)
    ])
def test_get_set(model, case_name, external_source_vectors, external_source_rate,
                 timesteps):
    """Tests the get/set methods"""

    op = CoupledOperator(model, CHAIN_PATH)
    number_of_timesteps = 2
    transfer = ExternalSourceRates(op, model, number_of_timesteps)

    if timesteps is None:
        timesteps = np.arange(number_of_timesteps)

    # Test by Openmc material, material name and material id
    material= [m for m in model.materials if m.depletable][0]

    for material_input in [material, material.name, material.id]:
        for external_source_vector in external_source_vectors:
            if case_name == 'rates_invalid_1':
                with pytest.raises(ValueError, match='Gb is not a valid '
                                   'nuclide or element.'):
                    transfer.set_external_source_rate(material_input,
                                               external_source_vector,
                                               external_source_rate)
            elif case_name == 'rates_invalid_2':
                with pytest.raises(ValueError, match='Cannot add element Pu '
                                             'as it is not naturally abundant. '
                                             'Specify a nuclide vector instead. '):
                    transfer.set_external_source_rate(material_input,
                                               external_source_vector,
                                               external_source_rate)
            else:
                transfer.set_external_source_rate(material_input,
                                               external_source_vector,
                                               external_source_rate,
                                               timesteps=timesteps)
                for component, percent in external_source_vector.items():
                    split_component = re.split(r'\d+', component)
                    if len(split_component) == 1:
                        for nuc,frac in openmc.data.isotopes(component):
                            assert transfer.get_external_rate(
                                material_input, nuc)[0] == pytest.approx(
                                                    external_source_rate \
                                                    * percent \
                                                    * frac \
                                                    * openmc.data.AVOGADRO \
                                                    / openmc.data.atomic_mass(nuc))
                    else:
                        assert transfer.get_external_rate(
                            material_input, component)[0] == pytest.approx(
                                                    external_source_rate \
                                                    * percent \
                                                    * openmc.data.AVOGADRO \
                                                    / openmc.data.atomic_mass(component))
                        assert component in transfer.get_components(material_input)
                assert np.all(transfer.get_material_timesteps(
                        material_input) == timesteps)


@pytest.mark.parametrize("external_source_rate_units, unit_conv", [
    ('g/s', 1),
    ('g/sec', 1),
    ('g/min', _SECONDS_PER_MINUTE),
    ('g/minute', _SECONDS_PER_MINUTE),
    ('g/h', _SECONDS_PER_HOUR),
    ('g/hr', _SECONDS_PER_HOUR),
    ('g/hour', _SECONDS_PER_HOUR),
    ('g/d', _SECONDS_PER_DAY),
    ('g/day', _SECONDS_PER_DAY),
    ('g/a', _SECONDS_PER_JULIAN_YEAR),
    ('g/year', _SECONDS_PER_JULIAN_YEAR),
    ])
def test_units(external_source_rate_units, unit_conv, model):
    """ Units testing"""
    # create external rate Xe
    components = ['Xe135', 'U235']
    external_source_rate = 1.0
    number_of_timesteps = 2
    op = CoupledOperator(model, CHAIN_PATH)
    transfer = ExternalSourceRates(op, model, number_of_timesteps)

    for component in components:
        transfer.set_external_source_rate('f', {component:1},
                                external_source_rate * unit_conv \
                                * openmc.data.atomic_mass(component) \
                                /openmc.data.AVOGADRO ,
                                external_source_rate_units=external_source_rate_units)
        assert transfer.get_external_rate('f', component)[0] == pytest.approx(
                                                        external_source_rate)


def test_external_source(run_in_tmpdir, model):
    """Tests external source depletion class without neither reaction rates nor
    decay but only external source rates"""
    # create transfer rate for U
    vector = {'U235':1}
    external_source = 10 #grams
    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_external_source_rate('f', vector, external_source/(24*3600))
    integrator.integrate()

    # Get number of U238 atoms from results
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(model.materials[0], "U235")

    # Ensure number of atoms equal external source
    assert atoms[1]-atoms[0] == pytest.approx(external_source \
                                             *openmc.data.AVOGADRO \
                                             /openmc.data.atomic_mass('U235'))
    assert atoms[2]-atoms[1] == pytest.approx(external_source \
                                             *openmc.data.AVOGADRO \
                                             /openmc.data.atomic_mass('U235'))
