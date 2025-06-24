"""Test one-group cross section generation"""
from pathlib import Path

import numpy as np
import pytest
import openmc
from openmc.deplete import MicroXS, get_microxs_and_flux

from tests.regression_tests import config

CHAIN_FILE = Path(__file__).parents[2] / "chain_simple.xml"

@pytest.fixture(scope="module")
def model():
    fuel = openmc.Material(name="uo2")
    fuel.add_nuclide("U235", 1.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cc", 10.4)

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=fuel)
    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.particles = 1000
    settings.inactive = 5
    settings.batches = 10

    return openmc.Model(geometry, settings=settings)


@pytest.mark.parametrize(
    "domain_type, rr_mode",
    [
        ("materials", "direct"),
        ("materials", "flux"),
        ("mesh", "direct"),
        ("mesh", "flux"),
    ]
)
def test_from_model(model, domain_type, rr_mode):
    if domain_type == 'materials':
        domains = list(model.geometry.get_all_materials().values())
    elif domain_type == 'mesh':
        mesh = openmc.RegularMesh()
        mesh.lower_left = (-10., -10.)
        mesh.upper_right = (10., 10.)
        mesh.dimension = (1, 1)
        domains = mesh
    nuclides = ['U235', 'O16', 'Xe135']
    kwargs = {
        'reaction_rate_mode': rr_mode,
        'chain_file': CHAIN_FILE,
        'path_statepoint': 'neutron_transport.h5',
    }
    if rr_mode == 'flux':
        kwargs['energies'] = 'CASMO-40'
    _, test_xs = get_microxs_and_flux(model, domains, nuclides, **kwargs)
    if config['update']:
        test_xs[0].to_csv(f'test_reference_{domain_type}_{rr_mode}.csv')

    # Make sure results match reference results
    ref_xs = MicroXS.from_csv(f'test_reference_{domain_type}_{rr_mode}.csv')
    np.testing.assert_allclose(test_xs[0].data, ref_xs.data, rtol=1e-11)

    # Make sure statepoint file was saved
    assert Path('neutron_transport.h5').exists()
