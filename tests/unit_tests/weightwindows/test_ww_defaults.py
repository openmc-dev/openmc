"""Default weight window energy bounds must follow the particle type."""

import numpy as np
import pytest
import openmc
import openmc.lib


@pytest.fixture
def lib_model(run_in_tmpdir):
    """Minimal model with both neutron and photon data available."""
    openmc.reset_auto_ids()
    model = openmc.Model()

    water = openmc.Material()
    water.set_density('g/cm3', 1.0)
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)

    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=water, region=-sphere)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = 'fixed source'
    model.settings.particles = 100
    model.settings.batches = 1
    model.settings.photon_transport = True

    model.export_to_model_xml()
    openmc.lib.init()
    # Required: initialize_data() runs here, not in openmc_init(), and it is
    # what narrows data::energy_min/max to the loaded data for each particle
    openmc.lib.simulation_init()
    yield model
    openmc.lib.simulation_finalize()
    openmc.lib.finalize()


def _lib_mesh():
    mesh = openmc.lib.RegularMesh()
    mesh.dimension = (2, 2, 2)
    mesh.set_parameters(lower_left=(-1.0, -1.0, -1.0),
                        upper_right=(1.0, 1.0, 1.0))
    return mesh


@pytest.mark.parametrize('particle', ('neutron', 'photon'))
def test_default_energy_bounds_follow_particle(lib_model, particle):
    """Defaults are derived after the particle type is known, not before."""
    ww = openmc.lib.WeightWindows(300 if particle == 'neutron' else 301)
    ww.mesh = _lib_mesh()
    ww.particle = particle

    bounds = np.asarray(ww.energy_bounds)
    assert bounds.size == 2

    # Compare against the range the library reports for this particle. Using
    # the other particle's range would be the symptom of deriving defaults in
    # the constructor.
    other = 'photon' if particle == 'neutron' else 'neutron'
    other_ww = openmc.lib.WeightWindows(400 if particle == 'neutron' else 401)
    other_ww.mesh = _lib_mesh()
    other_ww.particle = other
    other_bounds = np.asarray(other_ww.energy_bounds)

    assert not np.allclose(bounds, other_bounds), (
        f'{particle} and {other} weight windows have identical default energy '
        'bounds, which suggests the defaults were not derived from the '
        'particle type'
    )


def test_explicit_energy_bounds_survive_particle_type(lib_model):
    """Setting the particle type must not overwrite an explicit grid."""
    ww = openmc.lib.WeightWindows(302)
    ww.mesh = _lib_mesh()
    ww.energy_bounds = (1.0e3, 1.0e5, 1.0e7)
    ww.particle = 'photon'

    np.testing.assert_allclose(ww.energy_bounds, (1.0e3, 1.0e5, 1.0e7))


def test_bounds_survive_particle_type(lib_model):
    """Setting the particle type must not discard weight window bounds."""
    ww = openmc.lib.WeightWindows(303)
    ww.mesh = _lib_mesh()
    ww.energy_bounds = (0.0, 1.0e7)
    lower = np.arange(1.0, 9.0)
    ww.bounds = lower, 5.0 * lower
    ww.particle = 'photon'

    np.testing.assert_allclose(ww.bounds[0], lower)
    np.testing.assert_allclose(ww.bounds[1], 5.0 * lower)
