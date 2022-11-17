import numpy as np
import openmc
from tests.testing_harness import PyAPITestHarness

def compute_angular_distribution(cfg, E0, N):
    """Return a openmc.model.Model() object for a monoenergetic pencil
     beam hitting a 1 mm sphere filled with the material defined by
     the cfg string, and compute the angular distribution"""

    # Material definition

    m1 = openmc.Material.from_ncrystal(cfg)
    materials = openmc.Materials([m1])

    # Geometry definition

    sample_sphere = openmc.Sphere(r=0.1)
    outer_sphere = openmc.Sphere(r=100, boundary_type="vacuum")
    cell1 = openmc.Cell(region= -sample_sphere, fill=m1)
    cell2_region= +sample_sphere&-outer_sphere
    cell2 = openmc.Cell(region= cell2_region, fill=None)
    uni1 = openmc.Universe(cells=[cell1, cell2])
    geometry = openmc.Geometry(uni1)

    # Source definition

    source = openmc.Source()
    source.space = openmc.stats.Point(xyz = (0,0,-20))
    source.angle = openmc.stats.Monodirectional(reference_uvw = (0,0,1))
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
    filter1 = openmc.filter.SurfaceFilter(sample_sphere)
    filter2 = openmc.filter.PolarFilter(np.linspace(0,np.pi,180+1))
    filter3 = openmc.filter.CellFromFilter(cell1)
    tally1.filters = [filter1, filter2, filter3]
    tallies = openmc.Tallies([tally1])

    return openmc.model.Model(geometry, materials, settings, tallies)


class NCrystalTest(PyAPITestHarness):
    def _get_results(self):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)
        tal = sp.get_tally(name='angular distribution')
        df = tal.get_pandas_dataframe()
        return df.to_string()

def test_ncrystal():
    NParticles = 100000
    T = 293.6 # K
    E0 = 0.012 # eV
    cfg = 'Al_sg225.ncmat'
    test = compute_angular_distribution(cfg, E0, NParticles)
    harness = NCrystalTest('statepoint.10.h5', model=test)
    harness.main()
