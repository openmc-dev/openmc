import os
import shutil

import h5py
import numpy as np
import pytest
import openmc

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


@pytest.fixture
def model(request):
    openmc.reset_auto_ids()
    marker = request.node.get_closest_marker("surf_source_op")
    surf_source_op = marker.args[0]

    openmc_model = openmc.model.Model()

    # Materials
    # None

    # Geometry
    # Concentric void spheres
    # - Innermost sphere to bank surface sources
    # - Second shell to tally cell flux
    # - Outermost sphere as vacuum boundary
    sph_1 = openmc.Sphere(r=1.0)  # Surface to bank/write sources.
    sph_2 = openmc.Sphere(r=2.0)
    sph_3 = openmc.Sphere(r=2.5)
    sph_4 = openmc.Sphere(r=4.0, boundary_type="vacuum")
    cell_1 = openmc.Cell(region=-sph_1)
    cell_2 = openmc.Cell(region=+sph_1 & -sph_2)
    cell_3 = openmc.Cell(region=+sph_2 & -sph_3)  # Cell to tally flux.
    cell_4 = openmc.Cell(region=+sph_3 & -sph_4)
    root = openmc.Universe(cells=[cell_1, cell_2, cell_3, cell_4])
    openmc_model.geometry = openmc.Geometry(root)

    # Settings
    openmc_model.settings.run_mode = "fixed source"
    openmc_model.settings.particles = 1000
    openmc_model.settings.batches = 10
    openmc_model.settings.seed = 1

    if surf_source_op == "write":
        point = openmc.stats.Point((0, 0, 0))
        pt_src = openmc.IndependentSource(space=point)
        openmc_model.settings.source = pt_src

        openmc_model.settings.surf_source_write = {
            "surface_ids": [1],
            "max_particles": 1000,
        }
    elif surf_source_op == "read":
        openmc_model.settings.surf_source_read = {"path": "surface_source_true.h5"}

    # Tallies
    tal = openmc.Tally()
    cell_filter = openmc.CellFilter(cell_3)
    tal.filters = [cell_filter]
    tal.scores = ["flux"]
    openmc_model.tallies.append(tal)

    return openmc_model


class SurfaceSourceTestHarness(PyAPITestHarness):
    def _test_output_created(self):
        """Make sure surface_source.h5 has also been created."""
        super()._test_output_created()
        # Check if 'surface_source.h5' has been created.
        if self._model.settings.surf_source_write:
            assert os.path.exists(
                "surface_source.h5"
            ), "Surface source file does not exist."

    def _compare_output(self):
        """Make sure the current surface_source.h5 agree with the reference."""
        if self._model.settings.surf_source_write:
            with h5py.File("surface_source_true.h5", "r") as f:
                source_true = f["source_bank"][()]
                # Convert dtye from mixed to a float for comparison assertion
                source_true.dtype = "float64"
            with h5py.File("surface_source.h5", "r") as f:
                source_test = f["source_bank"][()]
                # Convert dtye from mixed to a float for comparison assertion
                source_test.dtype = "float64"
            np.testing.assert_allclose(
                np.sort(source_true), np.sort(source_test), atol=1e-07
            )

    def execute_test(self):
        """Build input XMLs, run OpenMC, check output and results."""
        try:
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

    def _overwrite_results(self):
        """Overwrite the results_true with the results_test."""
        shutil.copyfile("results_test.dat", "results_true.dat")
        if os.path.exists("surface_source.h5"):
            shutil.copyfile("surface_source.h5", "surface_source_true.h5")

    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        super()._cleanup()
        fs = "surface_source.h5"
        if os.path.exists(fs):
            os.remove(fs)


@pytest.mark.surf_source_op("write")
def test_surface_source_write(model, monkeypatch):
    # Test result is based on 1 MPI process
    monkeypatch.setitem(config, "mpi_np", "1")
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, "inputs_true_write.dat"
    )
    harness.main()


@pytest.mark.surf_source_op("read")
def test_surface_source_read(model):
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, "inputs_true_read.dat"
    )
    harness.main()
