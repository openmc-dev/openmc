"""Test the 'collision_track' setting.

Results
-------

All results are generated using only 1 MPI process.

All results are generated using 1 thread except for "test_consistency_low_realization_number".
This specific test verifies that when the number of realization (i.e., point being candidate
to be stored) is lower than the capacity, results are reproducible even with multiple
threads (i.e., there is no potential thread competition that would produce different
results in that case).

All results are generated using the history-based mode except for cases e01 to e03.

All results are visually verified using the '_visualize.py' script in the regression test folder.

OpenMC models
-------------

Four OpenMC models with CSG-only geometries are used to cover the transmission, vacuum,
reflective and periodic Boundary Conditions (BC):

- model_1: cylindrical core in 2 boxes (vacuum and transmission BC),

# Test cases for simulation parameters using CSG-only geometries
# ============================================================
# Each test case is defined by a combination of folder name, model name, and specific parameters.
# Below is a summary of the parameters used in the test cases:
#
# - max_collisions: Maximum number of particles to track in the simulation.
# - reactions: List of MT numbers (reaction types- 2 for scattering, 18 for fission, 101 for absorbtion).
# - cell_ids: IDs of specific cells in the model.
# - mat_ids: Material IDs for filtering particles.
# - nuclide_ids: Nuclide IDs for filtering particles.
# - univ_ids: Universe IDs for filtering particles.
# - E_threshold: Energy threshold for filtering particles (optional).
#
# The test cases are designed to validate the behavior of the simulation under various configurations.

*: BC stands for Boundary Conditions, T for Transmission, R for Reflective, and V for Vacuum.

An additional case, called 'case-a01', is used to check that the results are comparable when
the number of threads is set to 2 if the number of realization is lower than the capacity.


*: BC stands for Boundary Conditions, T for Transmission, and V for Vacuum.

Notes:

- The test cases list is non-exhaustive compared to the number of possible combinations.
  Test cases have been selected based on use and internal code logic.



TODO:

- Test with a lattice.

"""

import os
import shutil
from pathlib import Path

import h5py
import numpy as np
import openmc
import openmc.lib
import pytest

from tests.testing_harness import PyAPITestHarness
from tests.testing_harness import CollisionTrackTestHarness
from tests.regression_tests import config


@pytest.fixture(scope="function")
def two_threads(monkeypatch):
    """Set the number of OMP threads to 2 for the test."""
    monkeypatch.setenv("OMP_NUM_THREADS", "2")


@pytest.fixture(scope="function")
def single_process(monkeypatch):
    """Set the number of MPI process to 1 for the test."""
    monkeypatch.setitem(config, "mpi_np", "1")


@pytest.fixture(scope="module")
def model_1():
    """Cylindrical core contained in a first box which is contained in a larger box.
    A lower universe is used to describe the interior of the first box which
    contains the core and its surrounding space.

    """
    openmc.reset_auto_ids()
    model = openmc.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material(material_id=1)
    fuel.add_nuclide("U234", 0.0004524)
    fuel.add_nuclide("U235", 0.0506068)
    fuel.add_nuclide("U238", 0.9487090)
    fuel.add_nuclide("U236", 0.0002318)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material(material_id=11)
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", 1.0)

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
    core_lower_plane = openmc.ZPlane(-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region, cell_id=22)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(
        fill=water, region=outside_core_region, cell_id=33)

    # Universe
    inside_box1_universe = openmc.Universe(
        cells=[core, outside_core], universe_id=77)

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_rpp = openmc.model.RectangularParallelepiped(
        -box1_size / 2.0, box1_size / 2.0,
        -box1_size / 2.0, box1_size / 2.0,
        -box1_size / 2.0, box1_size / 2.0,
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=-box1_rpp, cell_id=5)

    # -----------------------------------------------------------------------------
    # Box 2
    # -----------------------------------------------------------------------------

    # Parameters
    box2_size = 8

    # Surfaces
    box2_rpp = openmc.model.RectangularParallelepiped(
        -box2_size / 2.0, box2_size / 2.0,
        -box2_size / 2.0, box2_size / 2.0,
        -box2_size / 2.0, box2_size / 2.0,
        boundary_type="vacuum"
    )

    # Cell
    box2 = openmc.Cell(fill=water, region=-box2_rpp & +box1_rpp, cell_id=8)

    # Register geometry
    model.geometry = openmc.Geometry([box1, box2])

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.particles = 100
    model.settings.batches = 5
    model.settings.inactive = 1
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
        space=distribution, constraints={'fissionable': True})

    return model


@pytest.mark.parametrize(
    "folder, model_name, parameter",
    [("case_1_Reactions", "model_1", {"max_collisions": 300, "reactions": ["(n,fission)", 101]}),
        ("case_2_Cell_ID", "model_1", {
         "max_collisions": 300, "cell_ids": [22]}),
        ("case_3_Material_ID", "model_1", {
         "max_collisions": 300, "material_ids": [1]}),
        ("case_4_Nuclide_ID", "model_1", {
         "max_collisions": 300, "nuclide_ids": ["O16", "U235"]}),
        ("case_5_Universe_ID", "model_1", {
         "max_collisions": 300, "cell_ids": [22], "universe_ids": [77]}),
        ("case_6_deposited_energy_threshold", "model_1", {
         "max_collisions": 300, "deposited_E_threshold": 5.5e5}),
        ("case_7_all_parameters_used_together", "model_1", {
            "max_collisions": 300,
            "reactions": ["elastic", 18, "(n,disappear)"],
            "material_ids": [1, 11],
            "universe_ids": [77],
            "nuclide_ids": ["U238", "U235", "H1", "U234"],
            "cell_ids": [22, 33],
            "deposited_E_threshold": 1e5})
     ],
)
def test_collision_track_several_cases(
    folder, model_name, parameter, request
):
    # Since for these tests the actual number of collisions recorded is < max_collisions,
    # we can run them with 1 or 2 threads, and in history or event mode.
    model = request.getfixturevalue(model_name)
    model.settings.collision_track = parameter
    harness = CollisionTrackTestHarness(
        "statepoint.5.h5", model=model, workdir=folder
    )
    harness.main()


@pytest.mark.skipif(config["event"] is True, reason="Results from history-based mode.")
def test_collision_track_2threads(model_1, two_threads, single_process):
    # This test checks that the `max_collisions` setting is honored:
    # no collisions beyond the specified limit should be recorded.
    #
    # For the result to be reproducible, the number of threads and
    # the transport mode (history vs. event) must remain fixed.
    assert os.environ["OMP_NUM_THREADS"] == "2"
    assert config["mpi_np"] == "1"
    model_1.settings.collision_track = {
        "max_collisions": 200
    }
    harness = CollisionTrackTestHarness(
        "statepoint.5.h5", model=model_1, workdir="case_8_2threads"
    )
    harness.main()
