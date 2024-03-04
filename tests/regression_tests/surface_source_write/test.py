"""Test the 'surface_source_write' setting.

To avoid the impact of thread competition during the source point storage, this test is limited
to a single thread via the fixture single_thread. Without this configuration, there is a chance
that source files can be different for the same inputs if the number of realization (i.e., point
being candidate to be stored) is higher than the capacity.

Results are generated using only 1 MPI process.

Three OpenMC models are used to cover the transmission, vacuum and reflective boundary conditions:

- model_1: complete model with a cylindrical core in 2 boxes,
- model_2: simplified model with a cylindrical core in 1 box (vacuum boundary conditions),
- model_3: simplified model with a cylindrical core in 1 box (reflective boundary conditions).

Results are visually verified using the '_visualize.py' script in the regression test folder.

Test cases:

========  =======  =========  =========================  =====  ===================================
Folder    Model    Surface    Cell                       BC*    Expected particles
========  =======  =========  =========================  =====  ===================================
case-1    model_1  No         No                         T+V    Particles crossing any surface in
                                                                the model
case-2    model_1  1          No                         T      Particles crossing this surface
                                                                only
case-3    model_1  Multiple   No                         T      Particles crossing the declared
                                                                surfaces
case-4    model_1  Multiple   cell (lower universe)      T      Particles crossing the declared
                                                                surfaces that come from or are
                                                                coming to the cell
case-5    model_1  Multiple   cell (root universe)       T      Particles crossing the declared
                                                                surfaces that come from or are
                                                                coming to the cell
case-6    model_1  No         cell (lower universe)      T      Particles crossing any surface that
                                                                come from or are coming to the cell
case-7    model_1  No         cell (root universe)       T      Particles crossing any surface that
                                                                come from or are coming to the cell
case-8    model_1  No         cellfrom (lower universe)  T      Particles crossing any surface that
                                                                come from the cell
case-9    model_1  No         cellto (lower universe)    T      Particles crossing any surface that
                                                                are coming to the cell
case-10   model_1  No         cellfrom (root universe)   T      Particles crossing any surface that
                                                                come from the cell
case-11   model_1  No         cellto (root universe)     T      Particles crossing any surface that
                                                                are coming to the cell
case-12   model_2  Multiple   No                         V      Particles crossing the declared
                                                                surfaces
case-13   model_2  Multiple   cell (root universe)       V      Particles crossing any surface that
                                                                come from or are coming to the cell
case-14   model_2  Multiple   cellfrom (root universe)   V      Particles crossing any surface that
                                                                are coming to the cell
case-15   model_2  Multiple   cellto (root universe)     V      None
case-16   model_3  Multiple   No                         R      Particles crossing the declared
                                                                surfaces
case-17   model_3  Multiple   cell (root universe)       R      None
case-18   model_3  Multiple   cellfrom (root universe)   R      None
case-19   model_3  Multiple   cellto (root universe)     R      None
========  =======  =========  =========================  =====  ===================================

*: BC stands for Boundary Conditions, T for Transmission, R for Reflective, and V for Vacuum.

An additional case, called 'case-20', is used to check that the results are comparable when
the number of threads is set to 2 if the number of realization is lower than the capacity.

Notes:

- The test cases list is non-exhaustive compared to the number of possible combinations.
  Test cases have been selected based on use and internal code logic.
- Cases 8 to 11 are testing that the feature still works even if the level of coordinates
  before and after crossing a surface is different,
- Cases that should return an error are tested in the 'test_exceptions' unit test
  from 'test_surf_source_write.py'.

TODO:

- Test with a lattice,
- Test with mesh,
- Test with periodic boundary conditions.

"""

import os
import shutil

import h5py
import numpy as np
import pytest
import openmc

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


@pytest.fixture(scope="function")
def single_thread(monkeypatch):
    """Set the number of OMP threads to 1 for the test."""
    monkeypatch.setenv("OMP_NUM_THREADS", "1")


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
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
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
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0)
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0)
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0)
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0)
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0)
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0)

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # -----------------------------------------------------------------------------
    # Box 2
    # -----------------------------------------------------------------------------

    # Parameters
    box2_size = 8

    # Surfaces
    box2_lower_plane = openmc.ZPlane(z0=-box2_size / 2.0, boundary_type="vacuum")
    box2_upper_plane = openmc.ZPlane(z0=box2_size / 2.0, boundary_type="vacuum")
    box2_left_plane = openmc.XPlane(x0=-box2_size / 2.0, boundary_type="vacuum")
    box2_right_plane = openmc.XPlane(x0=box2_size / 2.0, boundary_type="vacuum")
    box2_rear_plane = openmc.YPlane(y0=-box2_size / 2.0, boundary_type="vacuum")
    box2_front_plane = openmc.YPlane(y0=box2_size / 2.0, boundary_type="vacuum")

    # Region
    inside_box2 = (
        +box2_lower_plane
        & -box2_upper_plane
        & +box2_left_plane
        & -box2_right_plane
        & +box2_rear_plane
        & -box2_front_plane
    )
    outside_box1 = (
        -box1_lower_plane
        | +box1_upper_plane
        | -box1_left_plane
        | +box1_right_plane
        | -box1_rear_plane
        | +box1_front_plane
    )

    box2_region = inside_box2 & outside_box1

    # Cell
    box2 = openmc.Cell(fill=water, region=box2_region)

    # Root universe
    root = openmc.Universe(cells=[box1, box2])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
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
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


@pytest.fixture
def model_2():
    """Cylindrical core contained in a box.
    A lower universe is used to describe the interior of the box which
    contains the core and its surrounding space.

    The box is defined with vacuum boundary conditions.

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
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
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0, boundary_type="vacuum")
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0, boundary_type="vacuum")
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0, boundary_type="vacuum")
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0, boundary_type="vacuum")
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0, boundary_type="vacuum")
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0, boundary_type="vacuum")

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # Root universe
    root = openmc.Universe(cells=[box1])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
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
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


@pytest.fixture
def model_3():
    """Cylindrical core contained in a box.
    A lower universe is used to describe the interior of the box which
    contains the core and its surrounding space.

    The box is defined with reflective boundary conditions.

    """
    openmc.reset_auto_ids()
    model = openmc.model.Model()

    # =============================================================================
    # Materials
    # =============================================================================

    fuel = openmc.Material()
    fuel.add_element("U", 1.0, enrichment=5.0)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", 11.0)

    water = openmc.Material()
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
    core_lower_plane = openmc.ZPlane(z0=-core_height / 2.0)
    core_upper_plane = openmc.ZPlane(z0=core_height / 2.0)

    # Region
    core_region = -core_cylinder & +core_lower_plane & -core_upper_plane

    # Cells
    core = openmc.Cell(fill=fuel, region=core_region)
    outside_core_region = +core_cylinder | -core_lower_plane | +core_upper_plane
    outside_core = openmc.Cell(fill=water, region=outside_core_region)

    # Universe
    inside_box1_universe = openmc.Universe(cells=[core, outside_core])

    # -----------------------------------------------------------------------------
    # Box 1
    # -----------------------------------------------------------------------------

    # Parameters
    box1_size = 6.0

    # Surfaces
    box1_lower_plane = openmc.ZPlane(z0=-box1_size / 2.0, boundary_type="reflective")
    box1_upper_plane = openmc.ZPlane(z0=box1_size / 2.0, boundary_type="reflective")
    box1_left_plane = openmc.XPlane(x0=-box1_size / 2.0, boundary_type="reflective")
    box1_right_plane = openmc.XPlane(x0=box1_size / 2.0, boundary_type="reflective")
    box1_rear_plane = openmc.YPlane(y0=-box1_size / 2.0, boundary_type="reflective")
    box1_front_plane = openmc.YPlane(y0=box1_size / 2.0, boundary_type="reflective")

    # Region
    box1_region = (
        +box1_lower_plane
        & -box1_upper_plane
        & +box1_left_plane
        & -box1_right_plane
        & +box1_rear_plane
        & -box1_front_plane
    )

    # Cell
    box1 = openmc.Cell(fill=inside_box1_universe, region=box1_region)

    # Root universe
    root = openmc.Universe(cells=[box1])

    # Register geometry
    model.geometry = openmc.Geometry(root)

    # =============================================================================
    # Settings
    # =============================================================================

    model.settings = openmc.Settings()
    model.settings.run_mode = "eigenvalue"
    model.settings.particles = 1000
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
    distribution = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    model.settings.source = openmc.source.IndependentSource(space=distribution)

    return model


def return_surface_source_data(filepath):
    """Read a surface source file and return a sorted array composed
    of flatten arrays of source data for each surface source point.

    TODO:

    - use read_source_file from source.py instead. Or a dedicated function
      to produce sorted list of source points for a given file.

    Parameters
    ----------
    filepath : str
        Path to the surface source file

    Returns
    -------
    data : np.array
        Sorted array composed of flatten arrays of source data for
        each surface source point

    """
    data = []
    keys = []

    # Read source file
    with h5py.File(filepath, "r") as f:
        for point in f["source_bank"]:
            r = point["r"]
            u = point["u"]
            e = point["E"]
            time = point["time"]
            wgt = point["wgt"]
            delayed_group = point["delayed_group"]
            surf_id = point["surf_id"]
            particle = point["particle"]

            key = (
                f"{r[0]:.10e} {r[1]:.10e} {r[2]:.10e} {u[0]:.10e} {u[1]:.10e} {u[2]:.10e}"
                f"{e:.10e} {time:.10e} {wgt:.10e} {delayed_group} {surf_id} {particle}"
            )

            keys.append(key)

            values = [*r, *u, e, time, wgt, delayed_group, surf_id, particle]
            assert len(values) == 12
            data.append(values)

    data = np.array(data)
    keys = np.array(keys)
    sorted_idx = np.argsort(keys)

    data = data[sorted_idx]

    return data


class SurfaceSourceWriteTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, model=None, inputs_true=None, workdir=None):
        super().__init__(statepoint_name, model, inputs_true)
        self.workdir = workdir

    def _test_output_created(self):
        """Make sure surface_source.h5 has also been created."""
        super()._test_output_created()
        if self._model.settings.surf_source_write:
            assert os.path.exists(
                "surface_source.h5"
            ), "Surface source file has not been created."

    def _compare_output(self):
        """Compare surface_source.h5 files."""
        if self._model.settings.surf_source_write:
            source_true = return_surface_source_data("surface_source_true.h5")
            source_test = return_surface_source_data("surface_source.h5")
            np.testing.assert_allclose(source_true, source_test, rtol=1e-07)

    def main(self):
        """Accept commandline arguments and either run or update tests."""
        if config["build_inputs"]:
            self.build_inputs()
        elif config["update"]:
            self.update_results()
        else:
            self.execute_test()

    def build_inputs(self):
        """Build inputs."""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
        finally:
            os.chdir(base_dir)

    def execute_test(self):
        """Build inputs, run OpenMC, and verify correct results."""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._compare_inputs()
            self._run_openmc()
            self._test_output_created()
            self._compare_output()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()
            os.chdir(base_dir)

    def update_results(self):
        """Update results_true.dat and inputs_true.dat"""
        base_dir = os.getcwd()
        try:
            os.chdir(self.workdir)
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._overwrite_inputs()
            self._run_openmc()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._overwrite_results()
        finally:
            self._cleanup()
            os.chdir(base_dir)

    def _overwrite_results(self):
        """Also add the 'surface_source.h5' file during overwriting."""
        super()._overwrite_results()
        if os.path.exists("surface_source.h5"):
            shutil.copyfile("surface_source.h5", "surface_source_true.h5")

    def _cleanup(self):
        """Also remove the 'surface_source.h5' file while cleaning."""
        super()._cleanup()
        fs = "surface_source.h5"
        if os.path.exists(fs):
            os.remove(fs)


@pytest.mark.parametrize(
    "folder, parameter",
    [
        ("case-1", {"max_particles": 3000}),
        ("case-2", {"max_particles": 3000, "surface_ids": [4]}),
        ("case-3", {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]}),
        (
            "case-4",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 2},
        ),
        (
            "case-5",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        ("case-6", {"max_particles": 3000, "cell": 2}),
        ("case-7", {"max_particles": 3000, "cell": 3}),
        ("case-8", {"max_particles": 3000, "cellfrom": 2}),
        ("case-9", {"max_particles": 3000, "cellto": 2}),
        ("case-10", {"max_particles": 3000, "cellfrom": 3}),
        ("case-11", {"max_particles": 3000, "cellto": 3}),
    ],
)
def test_surface_source_cell_model_1(
    folder, parameter, model_1, single_thread, single_process
):
    """Test on a generic model with vacuum and transmission boundary conditions."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model_1.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model_1, workdir=folder
    )
    harness.main()


@pytest.mark.parametrize(
    "folder, parameter",
    [
        ("case-12", {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]}),
        (
            "case-13",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        (
            "case-14",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellfrom": 3},
        ),
        (
            "case-15",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellto": 3},
        ),
    ],
)
def test_surface_source_cell_model_2(
    folder, parameter, model_2, single_thread, single_process
):
    """Test specifically on vacuum surfaces."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model_2.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model_2, workdir=folder
    )
    harness.main()


@pytest.mark.parametrize(
    "folder, parameter",
    [
        ("case-16", {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9]}),
        (
            "case-17",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cell": 3},
        ),
        (
            "case-18",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellfrom": 3},
        ),
        (
            "case-19",
            {"max_particles": 3000, "surface_ids": [4, 5, 6, 7, 8, 9], "cellto": 3},
        ),
    ],
)
def test_surface_source_cell_model_3(
    folder, parameter, model_3, single_thread, single_process
):
    """Test specifically on reflective surfaces."""
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert config["mpi_np"] == "1"
    model_3.settings.surf_source_write = parameter
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model_3, workdir=folder
    )
    harness.main()


def test_consistency_low_realization_number(model_1, two_threads, single_process):
    """The objective is to test that the results produced, in a case where
    the number of potential realization (particle storage) is low
    compared to the capacity of storage, are still consistent.

    This configuration ensures that the competition between threads does not
    occur and that the content of the source file created can be compared.

    """
    assert os.environ["OMP_NUM_THREADS"] == "2"
    assert config["mpi_np"] == "1"
    model_1.settings.surf_source_write = {
        "max_particles": 200,
        "surface_ids": [2],
        "cellfrom": 2,
    }
    harness = SurfaceSourceWriteTestHarness(
        "statepoint.5.h5", model=model_1, workdir="case-20"
    )
    harness.main()