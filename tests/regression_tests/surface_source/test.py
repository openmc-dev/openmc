import os
import shutil

import h5py
import numpy as np
import pytest
import openmc

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


def mcpl_to_array(filepath):
    import mcpl

    source = []
    with mcpl.MCPLFile(filepath) as f:
        for p in f.particles:
            source.append(
                [
                    *tuple(p.position),
                    *tuple(p.direction),
                    1.0e6 * p.ekin,
                    1.0e-3 * p.time,
                    p.weight,
                    p.pdgcode,
                ]
            )
    return np.sort(np.array(source), axis=0)


def assert_structured_arrays_close(arr1, arr2, rtol=1e-5, atol=1e-8):
    assert arr1.dtype == arr2.dtype

    for field in arr1.dtype.names:
        data1, data2 = arr1[field], arr2[field]
        if data1.dtype.names:
            assert_structured_arrays_close(data1, data2, rtol=rtol, atol=atol)
        else:
            np.testing.assert_allclose(data1, data2, rtol=rtol, atol=atol)


@pytest.fixture
def model(request):
    openmc.reset_auto_ids()
    operation, file_format = request.node.get_closest_marker("params").args

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
    sph_4 = openmc.Sphere(r=4.0, boundary_type='vacuum')
    cell_1 = openmc.Cell(region=-sph_1)
    cell_2 = openmc.Cell(region=+sph_1&-sph_2)
    cell_3 = openmc.Cell(region=+sph_2&-sph_3)  # Cell to tally flux.
    cell_4 = openmc.Cell(region=+sph_3&-sph_4)
    root = openmc.Universe(cells=[cell_1, cell_2, cell_3, cell_4])
    openmc_model.geometry = openmc.Geometry(root)

    # Settings
    openmc_model.settings.run_mode = 'fixed source'
    openmc_model.settings.particles = 1000
    openmc_model.settings.batches = 10
    openmc_model.settings.seed = 1

    if operation == 'write':
        point = openmc.stats.Point((0, 0, 0))
        pt_src = openmc.IndependentSource(space=point)
        openmc_model.settings.source = pt_src

        surf_source_write_settings = {'surface_ids': [1],
                                      'max_particles': 1000}
        if file_format == "mcpl":
            surf_source_write_settings["mcpl"] = True

        openmc_model.settings.surf_source_write = surf_source_write_settings
    elif operation == 'read':
        openmc_model.settings.surf_source_read = {'path': f"surface_source_true.{file_format}"}

    # Tallies
    tal = openmc.Tally()
    cell_filter = openmc.CellFilter(cell_3)
    tal.filters = [cell_filter]
    tal.scores = ['flux']
    openmc_model.tallies.append(tal)

    return openmc_model


class SurfaceSourceTestHarness(PyAPITestHarness):
    def __init__(self, statepoint_name, model=None, inputs_true=None, file_format="h5"):
        super().__init__(statepoint_name, model, inputs_true)
        self.file_format = file_format

    def _test_output_created(self):
        """Make sure the surface_source file has also been created."""
        super()._test_output_created()
        if self._model.settings.surf_source_write:
            assert os.path.exists(f"surface_source.{self.file_format}"), \
                'Surface source file does not exist.'

    def _compare_output(self):
        """Make sure the current surface_source.h5 agree with the reference."""
        if self._model.settings.surf_source_write:
            if self.file_format == "h5":
                with h5py.File("surface_source_true.h5", 'r') as f:
                    source_true = np.sort(f['source_bank'][()])
                with h5py.File("surface_source.h5", 'r') as f:
                    source_test = np.sort(f['source_bank'][()])
                assert_structured_arrays_close(source_true, source_test, atol=1e-07)
            elif self.file_format == "mcpl":
                source_true = mcpl_to_array("surface_source_true.mcpl")
                source_test = mcpl_to_array("surface_source.mcpl")
                np.testing.assert_allclose(source_true, source_test, rtol=1e-5, atol=1e-7)

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
        shutil.copyfile('results_test.dat', 'results_true.dat')
        if os.path.exists(f"surface_source.{self.file_format}"):
            shutil.copyfile(f"surface_source.{self.file_format}", f"surface_source_true.{self.file_format}")

    def _cleanup(self):
        """Delete statepoints, tally, and test files."""
        super()._cleanup()
        fs = f"surface_source.{self.file_format}"
        if os.path.exists(fs):
            os.remove(fs)


@pytest.mark.params('write', 'h5')
def test_surface_source_write(model, monkeypatch, request):
    monkeypatch.setitem(config, "mpi_np", "1")  # Results generated with 1 MPI process
    operation, file_format = request.node.get_closest_marker("params").args
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, f"inputs_true_{operation}_{file_format}.dat",
        file_format=file_format
    )
    harness.main()


@pytest.mark.params('read', 'h5')
def test_surface_source_read(model, request):
    operation, file_format = request.node.get_closest_marker("params").args
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, f"inputs_true_{operation}_{file_format}.dat",
        file_format=file_format
    )
    harness.main()


@pytest.mark.skipif(shutil.which("mcpl-config") is None, reason="MCPL is not available.")
@pytest.mark.params('write', 'mcpl')
def test_surface_source_write_mcpl(model, monkeypatch, request):
    monkeypatch.setitem(config, "mpi_np", "1")  # Results generated with 1 MPI process
    operation, file_format = request.node.get_closest_marker("params").args
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, f"inputs_true_{operation}_{file_format}.dat",
        file_format=file_format
    )
    harness.main()


@pytest.mark.skipif(shutil.which("mcpl-config") is None, reason="MCPL is not available.")
@pytest.mark.params('read', 'mcpl')
def test_surface_source_read_mcpl(model, request):
    operation, file_format = request.node.get_closest_marker("params").args
    harness = SurfaceSourceTestHarness(
        "statepoint.10.h5", model, f"inputs_true_{operation}_{file_format}.dat",
        file_format=file_format
    )
    harness.main()
