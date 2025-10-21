from tests.testing_harness import TestHarness

import openmc

import pytest

pytestmark = pytest.mark.skipif(
    not openmc.lib._dagmc_enabled(),
    reason="DAGMC CAD geometry is not enabled.")


class DagmcWeightWindowsTestHarness(TestHarness):

    def __init__(self):
        super().__init__('')

    def execute_test(self):
        # this test merely needs to run OpenMC without error
        try:
            self._run_openmc()
        finally:
            self._cleanup()


def test_dagmc_weight_windows_near_boundary():
    """A regression test that ensures particle splitting near a boundary does
    not result in lost particles due to a stale DAGMC history on the particle
    object"""

    # This DAGMC model consists of three nested cubes. The innermost cube
    # contains a fusion neutron source. The two outer cubes are filled with
    # tungsten. The entire model has weight windows defined on a mesh such that
    # particles will be split as they move outward from the source. The
    # outermost cubes are very similary in size, the outer cube is just slightly
    # larger than the inner cube. This means that particles moving outward will
    # frequently cross the boundary between the two cubes right after being
    # split. This test ensures that no particles are lost due to a stale DAGMC
    # history after splitting.
    harness = DagmcWeightWindowsTestHarness()
    harness.main()
