from pathlib import Path
import subprocess
import sys

import openmc

from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness


class EmissionSpectrumHarness(PyAPITestHarness):
    def _run_openmc(self):
        command = [sys.executable, str(Path(__file__).with_name('worker.py'))]
        if config['event']:
            command.append('--event')
        if config['mpi']:
            command = [config['mpiexec'], '-n', str(config['mpi_np'])] + command
        subprocess.run(command, check=True)


def test_runtime_emission_spectrum():
    model = openmc.examples.pwr_pin_cell()
    model.settings.particles = 1000
    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.generations_per_batch = 1
    model.settings.seed = 19
    model.settings.verbosity = 1
    tally = openmc.Tally()
    tally.filters = [openmc.EnergyFilter([0.0, 1.0, 1.e5, 1.e6, 20.e6])]
    tally.scores = ['flux', 'fission']
    model.tallies.append(tally)
    EmissionSpectrumHarness('statepoint.4.h5', model).main()
