import os
import openmc

from openmc.examples import random_ray_three_region_cube_with_detectors
from openmc.utility_funcs import change_directory
import pytest

from tests.regression_tests import config
from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def make_model():
    model = random_ray_three_region_cube_with_detectors()

    detector1_cells = model.geometry.get_cells_by_name("detector 1")
    detector2_cells = model.geometry.get_cells_by_name("detector 2")
    detector_cells = detector1_cells + detector2_cells

    strengths = [1.0]
    midpoints = [100.0]
    energy_distribution = openmc.stats.Discrete(x=midpoints, p=strengths)

    adj_source = openmc.IndependentSource(energy=energy_distribution, constraints={
                                          'domains': detector_cells}, strength=3.14)
    
    model.settings.random_ray['adjoint'] = True
    model.settings.random_ray['adjoint_source'] = adj_source
    model.settings.random_ray['volume_estimator'] = 'naive'
    return model


def test_random_ray_adjoint_local():
    harness = MGXSTestHarness('statepoint.10.h5', make_model())
    harness.main()


def test_hand_off_past_lattice_plane(tmp_path):
    """At 3 ranks, rays reach a rank boundary just past a lattice plane, where
    one handed over would be lost. The partition depends on timing, so this
    raises a flag in most runs (6 of 8 with the check removed), not in all."""
    if not config['mpi']:
        pytest.skip("Needs MPI.")
    if config['build_inputs']:
        return

    with change_directory(tmp_path):
        model = make_model()
    model.run(cwd=tmp_path, openmc_exec=config['exe'],
              mpi_args=[config['mpiexec'], '-n', '3'])
