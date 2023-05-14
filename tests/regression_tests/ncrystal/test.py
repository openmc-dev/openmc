from math import pi

import numpy as np
import openmc
import openmc.lib
import pytest

from tests.testing_harness import PyAPITestHarness

pytestmark = pytest.mark.skipif(
    not openmc.lib._ncrystal_enabled(),
    reason="NCrystal materials are not enabled.")


def pencil_beam_model(cfg, E0, N):
    """Return an openmc.Model() object for a monoenergetic pencil
     beam hitting a 1 mm sphere filled with the material defined by
     the cfg string, and compute the angular distribution"""

    # Material definition

    m1 = openmc.Material.from_ncrystal(cfg)
    materials = openmc.Materials([m1])

    # Geometry definition

    sample_sphere = openmc.Sphere(r=0.1)
    outer_sphere = openmc.Sphere(r=100, boundary_type="vacuum")
    cell1 = openmc.Cell(region=-sample_sphere, fill=m1)
    cell2_region = +sample_sphere & -outer_sphere
    cell2 = openmc.Cell(region=cell2_region, fill=None)
    geometry = openmc.Geometry([cell1, cell2])

    # Source definition

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, -20))
    source.angle = openmc.stats.Monodirectional(reference_uvw=(0, 0, 1))
    source.energy = openmc.stats.Discrete([E0], [1.0])

    # Execution settings

    settings = openmc.Settings()
    settings.source = source
    settings.run_mode = "fixed source"
    settings.batches = 10
    settings.particles = N

    # Tally definition

    tally1 = openmc.Tally(name="angular distribution")
    tally1.scores = ["current"]
    filter1 = openmc.SurfaceFilter(sample_sphere)
    filter2 = openmc.PolarFilter(np.linspace(0, pi, 180+1))
    filter3 = openmc.CellFromFilter(cell1)
    tally1.filters = [filter1, filter2, filter3]
    tallies = openmc.Tallies([tally1])

    return openmc.Model(geometry, materials, settings, tallies)


class NCrystalTest(PyAPITestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        with openmc.StatePoint(self._sp_name) as sp:
            tal = sp.get_tally(name='angular distribution')
            df = tal.get_pandas_dataframe()
        return df.to_string()


def test_ncrystal():
    n_particles = 100000
    T = 293.6  # K
    E0 = 0.012  # eV
    cfg = 'Al_sg225.ncmat'
    test = pencil_beam_model(cfg, E0, n_particles)
    harness = NCrystalTest('statepoint.10.h5', model=test)
    harness.main()
