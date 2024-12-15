"""This test ensures that the CellFromFilter works correctly even if the level of
coordinates (number of encapsulated universes) is different in the cell from
where the particle originates compared to the cell where the particle is going.

A matrix of reaction rates based on where the particle is coming from and
where it goes to is calculated and compared to the total reaction rate of the problem.
The components of this matrix are also compared to other components using symmetric
properties.

TODO:

- Test with a lattice,
- Test with mesh,
- Test with reflective boundary conditions,
- Test with periodic boundary conditions.

"""

from numpy.testing import assert_allclose, assert_equal
import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


RTOL = 1.0e-7
ATOL = 0.0


@pytest.fixture
def model():
    """Cylindrical core contained in a first box which is contained in a larger box.
    A lower universe is used to describe the interior of the first box which
    contains the core and its surrounding space."""
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_nuclide("U234", 0.0004524)
    fuel.add_nuclide("U235", 0.0506068)
    fuel.add_nuclide("U238", 0.9487090)
    fuel.add_nuclide("U236", 0.0002318)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 10.97)

    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)

    air = openmc.Material()
    air.add_element("O", 0.2)
    air.add_element("N", 0.8)
    air.set_density("g/cm3", 0.001225)

    # =============================================================================
    # Geometry
    # =============================================================================

    # -----------------------------------------------------------------------------
    # Cylindrical core
    # -----------------------------------------------------------------------------

    # Parameters
    core_radius = 2.0
    core_height = 4.0

    # Surfaces
    core_cylinder = openmc.ZCylinder(r=core_radius)
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=air, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 4.1

    # Surfaces
    box1_rpp = openmc.model.RectangularParallelepiped(
        -box1_size / 2.0,
        box1_size / 2.0,
        -box1_size / 2.0,
        box1_size / 2.0,
        -box1_size / 2.0,
        box1_size / 2.0,
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=-box1_rpp)

    # -----------------------------------------------------------------------------
    # Box 2
    # -----------------------------------------------------------------------------

    # Parameters
    box2_size = 12

    # Surfaces
    box2_rpp = openmc.model.RectangularParallelepiped(
        -box2_size / 2.0,
        box2_size / 2.0,
        -box2_size / 2.0,
        box2_size / 2.0,
        -box2_size / 2.0,
        box2_size / 2.0,
        boundary_type="vacuum",
    )

    # Cell
    box2 = openmc.Cell(fill=water, region=-box2_rpp & +box1_rpp)

    # Register geometry
    model.geometry = openmc.Geometry([box1, box2])

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.particles = 2000
    model.settings.batches = 15
    model.settings.inactive = 5
    model.settings.seed = 1

    bounds = [
        -core_radius,
        -core_radius,
        -core_height / 2.0,
        core_radius,
        core_radius,
        core_height / 2.0,
    ]
    distribution = openmc.stats.Box(bounds[:3], bounds[3:])
    model.settings.source = openmc.IndependentSource(
        space=distribution, constraints={"fissionable": True}
    )

    # =============================================================================
    # Tallies
    # =============================================================================

    in_core_filter = openmc.CellFilter([core])
    in_outside_core_filter = openmc.CellFilter([outside_core])
    in_box1_filter = openmc.CellFilter([box1])
    in_box2_filter = openmc.CellFilter([box2])

    from_core_filter = openmc.CellFromFilter([core])
    from_outside_core_filter = openmc.CellFromFilter([outside_core])
    from_box1_filter = openmc.CellFromFilter([box1])
    from_box2_filter = openmc.CellFromFilter([box2])

    t1_1 = openmc.Tally(name="total from 1 in 1")
    t1_1.filters = [from_core_filter, in_core_filter]
    t1_1.scores = ["total"]

    t1_2 = openmc.Tally(name="total from 1 in 2")
    t1_2.filters = [from_core_filter, in_outside_core_filter]
    t1_2.scores = ["total"]

    t1_3 = openmc.Tally(name="total from 1 in 3")
    t1_3.filters = [from_core_filter, in_box1_filter]
    t1_3.scores = ["total"]

    t1_4 = openmc.Tally(name="total from 1 in 4")
    t1_4.filters = [from_core_filter, in_box2_filter]
    t1_4.scores = ["total"]

    t2_1 = openmc.Tally(name="total from 2 in 1")
    t2_1.filters = [from_outside_core_filter, in_core_filter]
    t2_1.scores = ["total"]

    t2_2 = openmc.Tally(name="total from 2 in 2")
    t2_2.filters = [from_outside_core_filter, in_outside_core_filter]
    t2_2.scores = ["total"]

    t2_3 = openmc.Tally(name="total from 2 in 3")
    t2_3.filters = [from_outside_core_filter, in_box1_filter]
    t2_3.scores = ["total"]

    t2_4 = openmc.Tally(name="total from 2 in 4")
    t2_4.filters = [from_outside_core_filter, in_box2_filter]
    t2_4.scores = ["total"]

    t3_1 = openmc.Tally(name="total from 3 in 1")
    t3_1.filters = [from_box1_filter, in_core_filter]
    t3_1.scores = ["total"]

    t3_2 = openmc.Tally(name="total from 3 in 2")
    t3_2.filters = [from_box1_filter, in_outside_core_filter]
    t3_2.scores = ["total"]

    t3_3 = openmc.Tally(name="total from 3 in 3")
    t3_3.filters = [from_box1_filter, in_box1_filter]
    t3_3.scores = ["total"]

    t3_4 = openmc.Tally(name="total from 3 in 4")
    t3_4.filters = [from_box1_filter, in_box2_filter]
    t3_4.scores = ["total"]

    t4_1 = openmc.Tally(name="total from 4 in 1")
    t4_1.filters = [from_box2_filter, in_core_filter]
    t4_1.scores = ["total"]

    t4_2 = openmc.Tally(name="total from 4 in 2")
    t4_2.filters = [from_box2_filter, in_outside_core_filter]
    t4_2.scores = ["total"]

    t4_3 = openmc.Tally(name="total from 4 in 3")
    t4_3.filters = [from_box2_filter, in_box1_filter]
    t4_3.scores = ["total"]

    t4_4 = openmc.Tally(name="total from 4 in 4")
    t4_4.filters = [from_box2_filter, in_box2_filter]
    t4_4.scores = ["total"]

    tglobal = openmc.Tally(name="total")
    tglobal.scores = ["total"]

    model.tallies += [
        t1_1,
        t1_2,
        t1_3,
        t1_4,
        t2_1,
        t2_2,
        t2_3,
        t2_4,
        t3_1,
        t3_2,
        t3_3,
        t3_4,
        t4_1,
        t4_2,
        t4_3,
        t4_4,
        tglobal,
    ]
    return model


class CellFromFilterTest(PyAPITestHarness):

    def _compare_results(self):
        """Additional unit tests on the tally results to check
        consistency of CellFromFilter."""
        with openmc.StatePoint(self.statepoint_name) as sp:

            t1_1 = sp.get_tally(name="total from 1 in 1").mean
            t1_2 = sp.get_tally(name="total from 1 in 2").mean
            t1_3 = sp.get_tally(name="total from 1 in 3").mean
            t1_4 = sp.get_tally(name="total from 1 in 4").mean

            t2_1 = sp.get_tally(name="total from 2 in 1").mean
            t2_2 = sp.get_tally(name="total from 2 in 2").mean
            t2_3 = sp.get_tally(name="total from 2 in 3").mean
            t2_4 = sp.get_tally(name="total from 2 in 4").mean

            t3_1 = sp.get_tally(name="total from 3 in 1").mean
            t3_2 = sp.get_tally(name="total from 3 in 2").mean
            t3_3 = sp.get_tally(name="total from 3 in 3").mean
            t3_4 = sp.get_tally(name="total from 3 in 4").mean

            t4_1 = sp.get_tally(name="total from 4 in 1").mean
            t4_2 = sp.get_tally(name="total from 4 in 2").mean
            t4_3 = sp.get_tally(name="total from 4 in 3").mean
            t4_4 = sp.get_tally(name="total from 4 in 4").mean

            tglobal = sp.get_tally(name="total").mean

            # From 1 and 2 is equivalent to from 3
            assert_allclose(t1_1 + t2_1, t3_1, rtol=RTOL, atol=ATOL)
            assert_allclose(t1_2 + t2_2, t3_2, rtol=RTOL, atol=ATOL)
            assert_allclose(t1_3 + t2_3, t3_3, rtol=RTOL, atol=ATOL)
            assert_allclose(t1_4 + t2_4, t3_4, rtol=RTOL, atol=ATOL)

            # In 1 and 2 equivalent to in 3
            assert_allclose(t1_1 + t1_2, t1_3, rtol=RTOL, atol=ATOL)
            assert_allclose(t2_1 + t2_2, t2_3, rtol=RTOL, atol=ATOL)
            assert_allclose(t3_1 + t3_2, t3_3, rtol=RTOL, atol=ATOL)
            assert_allclose(t4_1 + t4_2, t4_3, rtol=RTOL, atol=ATOL)

            # Comparison to global from 3
            assert_allclose(t3_3 + t3_4 + t4_3 + t4_4, tglobal, rtol=RTOL, atol=ATOL)

            # Comparison to global from 1 and 2
            t_from_1_wo_3 = t1_1 + t1_2 + t1_4
            t_from_2_wo_3 = t2_1 + t2_2 + t2_4
            t_from_4_wo_3 = t4_1 + t4_2 + t4_4
            assert_allclose(
                t_from_1_wo_3 + t_from_2_wo_3 + t_from_4_wo_3,
                tglobal,
                rtol=RTOL,
                atol=ATOL,
            )

            # 1 cannot contribute to 4 and 4 cannot contribute to 1 by symmetry
            assert_equal(t1_4, np.zeros_like(t1_4))
            assert_equal(t4_1, np.zeros_like(t4_1))

        return super()._compare_results()


def test_filter_cellfrom(model):
    harness = CellFromFilterTest("statepoint.15.h5", model)
    harness.main()
