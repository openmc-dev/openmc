""" Tests for TransferRates class """

from pathlib import Path
from math import exp

import pytest
import numpy as np

import openmc
from openmc.deplete import CoupledOperator
from openmc.deplete.transfer_rates import TransferRates
from openmc.deplete.abc import (_SECONDS_PER_MINUTE, _SECONDS_PER_HOUR,
                                _SECONDS_PER_DAY, _SECONDS_PER_JULIAN_YEAR)

CHAIN_PATH = Path(__file__).parents[1] / "chain_simple.xml"

@pytest.fixture
def model():
    openmc.reset_auto_ids()
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
    h.depletable = True

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

@pytest.mark.parametrize("case_name, transfer_rates, timesteps", [
    ('elements', {'U': 0.01, 'Xe': 0.1}, None),
    ('nuclides', {'I135': 0.01, 'Gd156': 0.1, 'Gd157': 0.01}, None),
    ('nuclides_elements', {'I135': 0.01, 'Gd156': 0.1, 'Gd157': 0.01, 'U': 0.01,
                           'Xe': 0.1}, None),
    ('elements_nuclides', {'U': 0.01, 'Xe': 0.1, 'I135': 0.01, 'Gd156': 0.1,
                           'Gd157': 0.01}, None),
    ('multiple_transfer', {'U': 0.01, 'Xe': 0.1, 'I135': 0.01, 'Gd156': 0.1,
                           'Gd157': 0.01}, None),
    ('timesteps', {'U': 0.01, 'Xe': 0.1}, [1]),
    ('rates_invalid_1', {'Gd': 0.01, 'Gd157': 0.01, 'Gd156': 0.01}, None),
    ('rates_invalid_2', {'Gd156': 0.01, 'Gd157': 0.01, 'Gd': 0.01}, None),
    ('rates_invalid_3', {'Gb156': 0.01}, None),
    ('rates_invalid_4', {'Gb': 0.01}, None)
    ])
def test_get_set(model, case_name, transfer_rates, timesteps):
    """Tests the get/set methods"""
    op = CoupledOperator(model, CHAIN_PATH)
    number_of_timesteps = 2
    transfer = TransferRates(op, model.materials, number_of_timesteps)

    if timesteps is None:
        timesteps = np.arange(number_of_timesteps)

    # Test by Openmc material, material name and material id
    material, dest_material, dest_material2 = [m for m in model.materials
                                               if m.depletable]
    for material_input in [material, material.name, material.id]:
        for dest_material_input in [None, dest_material, dest_material.name,
                                    dest_material.id]:
            if case_name == 'rates_invalid_1':
                with pytest.raises(ValueError, match='Cannot add transfer '
                                    'rate for nuclide Gd157 to material 1 '
                                    'where element Gd already has a '
                                    'transfer rate.'):
                    for component, transfer_rate in transfer_rates.items():
                        transfer.set_transfer_rate(material_input,
                                                   [component],
                                                   transfer_rate)
            elif case_name == 'rates_invalid_2':
                with pytest.raises(ValueError, match='Cannot add transfer '
                                   f'rate for element Gd to material 1 with '
                                   r'transfer rate\(s\) for nuclide\(s\) '
                                   'Gd156, Gd157.'):
                    for component, transfer_rate in transfer_rates.items():
                        transfer.set_transfer_rate(material_input,
                                                   [component],
                                                   transfer_rate)
            elif case_name == 'rates_invalid_3':
                with pytest.raises(ValueError, match='Gb156 is not a valid '
                                   'nuclide or element.'):
                    for component, transfer_rate in transfer_rates.items():
                        transfer.set_transfer_rate(material_input,
                                                   [component],
                                                   transfer_rate)
            elif case_name == 'rates_invalid_4':
                with pytest.raises(ValueError, match='Gb is not a valid '
                                   'nuclide or element.'):
                    for component, transfer_rate in transfer_rates.items():
                        transfer.set_transfer_rate(material_input,
                                                   [component],
                                                   transfer_rate)
            else:
                for component, transfer_rate in transfer_rates.items():
                    transfer.set_transfer_rate(material_input, [component],
                                               transfer_rate,
                                               timesteps=timesteps,
                                               destination_material=\
                                               dest_material_input)
                    assert transfer.get_external_rate(
                        material_input, component, timesteps,
                        dest_material_input)[0] == transfer_rate
                    assert np.all(transfer.external_timesteps == timesteps)

                if timesteps is not None:
                    for timestep in timesteps:
                        assert transfer.get_components(material_input, timestep,
                            dest_material_input) == list(transfer_rates.keys())
                else:
                    assert transfer.get_components(material_input, timesteps,
                            dest_material_input) == list(transfer_rates.keys())

            if case_name == 'multiple_transfer':
                for dest2_material_input in [dest_material2, dest_material2.name,
                                    dest_material2.id]:
                    for component, transfer_rate in transfer_rates.items():
                        transfer.set_transfer_rate(material_input, [component],
                                               transfer_rate,
                                               destination_material=\
                                               dest2_material_input)
                        for id, dest_mat in zip([0,1],[dest_material,dest_material2]):
                            assert transfer.get_external_rate(
                                material_input, component, timesteps)[0] == transfer_rate

@pytest.mark.parametrize("transfer_rate_units, unit_conv", [
    ('1/s', 1),
    ('1/sec', 1),
    ('1/min', _SECONDS_PER_MINUTE),
    ('1/minute', _SECONDS_PER_MINUTE),
    ('1/h', _SECONDS_PER_HOUR),
    ('1/hr', _SECONDS_PER_HOUR),
    ('1/hour', _SECONDS_PER_HOUR),
    ('1/d', _SECONDS_PER_DAY),
    ('1/day', _SECONDS_PER_DAY),
    ('1/a', _SECONDS_PER_JULIAN_YEAR),
    ('1/year', _SECONDS_PER_JULIAN_YEAR),
    ])
def test_units(transfer_rate_units, unit_conv, model):
    """ Units testing"""
    # create transfer rate Xe
    components = ['Xe', 'U235']
    transfer_rate = 1e-5
    number_of_timesteps = 2
    op = CoupledOperator(model, CHAIN_PATH)
    transfer = TransferRates(op, model.materials, number_of_timesteps)

    for component in components:
        transfer.set_transfer_rate('f', [component], transfer_rate * unit_conv,
                                   transfer_rate_units=transfer_rate_units)
        for timestep in range(transfer.number_of_timesteps):
            assert transfer.get_external_rate('f', component, timestep)[0] == transfer_rate


def test_transfer(run_in_tmpdir, model):
    """Tests transfer depletion class without neither reaction rates nor decay
    but only transfer rates"""

    # create transfer rate for U
    element = ['U']
    transfer_rate = 1e-5
    op = CoupledOperator(model, CHAIN_PATH)
    integrator = openmc.deplete.PredictorIntegrator(
        op, [1,1], 0.0, timestep_units = 'd')
    integrator.add_transfer_rate('f', element, transfer_rate)
    integrator.integrate()

    # Get number of U238 atoms from results
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(model.materials[0], "U238")

    # Ensure number of atoms equal transfer decay
    assert atoms[1] == pytest.approx(atoms[0]*exp(-transfer_rate*3600*24))
    assert atoms[2] == pytest.approx(atoms[1]*exp(-transfer_rate*3600*24))

@pytest.mark.parametrize("case_name, buffer, ox", [
    ('redox', {'Gd157':1}, {'Gd': 3, 'U': 4}),
    ('buffer_invalid', {'Gd158':1}, {'Gd': 3, 'U': 4}),
    ('elm_invalid', {'Gd157':1}, {'Gb': 3, 'U': 4}),
    ])
def test_redox(case_name, buffer, ox, model):
    op = CoupledOperator(model, CHAIN_PATH)
    transfer = TransferRates(op, model)

    # Test by Openmc material, material name and material id
    material, dest_material, dest_material2 = [m for m in model.materials
                                               if m.depletable]
    for material_input in [material, material.name, material.id]:
        for dest_material_input in [dest_material, dest_material.name,
                                    dest_material.id]:

            if case_name == 'buffer_invalid':
                with pytest.raises(ValueError, match='Gd158 is not a valid '
                                   'nuclide.'):
                    transfer.set_redox(material_input, buffer, ox)

            elif case_name == 'elm_invalid':
                with pytest.raises(ValueError, match='Gb is not a valid '
                                   'element.'):
                    transfer.set_redox(material_input, buffer, ox)
            else:
                transfer.set_redox(material_input, buffer, ox)
                mat_id = transfer._get_material_id(material_input)
                assert transfer.redox[mat_id][0] == buffer
                assert transfer.redox[mat_id][1] == ox
