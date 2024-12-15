"""Test the meshborn filter using a fixed source calculation on a H1 sphere.

"""

from numpy.testing import assert_allclose
import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


RTOL = 1.0e-7
ATOL = 0.0


@pytest.fixture
def model():
    """Sphere of H1 with one hemisphere containing the source (x>0) and one
    hemisphere with no source (x<0).

    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Materials
    h1 = openmc.Material()
    h1.add_nuclide("H1", 1.0)
    h1.set_density("g/cm3", 1.0)
    model.materials = openmc.Materials([h1])

    # Core geometry
    r = 10.0
    sphere = openmc.Sphere(r=r, boundary_type="reflective")
    core = openmc.Cell(fill=h1, region=-sphere)
    model.geometry = openmc.Geometry([core])

    # Settings
    model.settings.run_mode = "fixed source"
    model.settings.particles = 2000
    model.settings.batches = 8
    distribution = openmc.stats.Box((0.0, -r, -r), (r, r, r))
    model.settings.source = openmc.IndependentSource(space=distribution)

    # Tallies
    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 2, 1)
    mesh.lower_left = (-r, -r, -r)
    mesh.upper_right = (r, r, r)

    f_1 = openmc.MeshFilter(mesh)
    f_2 = openmc.MeshBornFilter(mesh)

    t_1 = openmc.Tally(name="scatter")
    t_1.filters = [f_1, f_2]
    t_1.scores = ["scatter"]

    t_2 = openmc.Tally(name="scatter-mesh")
    t_2.filters = [f_1]
    t_2.scores = ["scatter"]

    t_3 = openmc.Tally(name="scatter-meshborn")
    t_3.filters = [f_2]
    t_3.scores = ["scatter"]

    t_4 = openmc.Tally(name="scatter-total")
    t_4.scores = ["scatter"]

    model.tallies = [t_1, t_2, t_3, t_4]

    return model


class MeshBornFilterTest(PyAPITestHarness):

    def _compare_results(self):
        """Additional unit tests on the tally results to check consistency."""
        with openmc.StatePoint(self.statepoint_name) as sp:

            t1 = sp.get_tally(name="scatter").mean.reshape(4, 4)
            t2 = sp.get_tally(name="scatter-mesh").mean.reshape(4)
            t3 = sp.get_tally(name="scatter-meshborn").mean.reshape(4)
            t4 = sp.get_tally(name="scatter-total").mean.reshape(1)

            # Consistency between mesh+meshborn matrix tally and meshborn tally
            for i in range(4):
                assert_allclose(t1[:, i].sum(), t3[i], rtol=RTOL, atol=ATOL)

            # Consistency between mesh+meshborn matrix tally and mesh tally
            for i in range(4):
                assert_allclose(t1[i, :].sum(), t2[i], rtol=RTOL, atol=ATOL)

            # Mesh cells in x<0 do not contribute to meshborn
            assert_allclose(t1[:, 0].sum(), np.zeros(4), rtol=RTOL, atol=ATOL)
            assert_allclose(t1[:, 2].sum(), np.zeros(4), rtol=RTOL, atol=ATOL)

            # Consistency with total scattering
            assert_allclose(t1.sum(), t4, rtol=RTOL, atol=ATOL)
            assert_allclose(t2.sum(), t4, rtol=RTOL, atol=ATOL)
            assert_allclose(t3.sum(), t4, rtol=RTOL, atol=ATOL)

        super()._compare_results()


def test_filter_meshborn(model):
    harness = MeshBornFilterTest("statepoint.8.h5", model)
    harness.main()
