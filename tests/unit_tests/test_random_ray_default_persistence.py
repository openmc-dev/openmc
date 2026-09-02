"""The random ray volume-estimator default must survive an in-process
finalize/re-initialize cycle (openmc.lib workflows such as iterative weight
window generation). The default is "auto", which is resolved to a concrete
estimator at the start of each run by overwriting the stored setting, so
openmc_finalize_random_ray() must restore "auto" between runs; if it leaves
the resolved value behind, later runs of the process inherit the previous
run's estimator instead of re-resolving. An adjoint run resolves to the
strict adaptive estimator, so running adjoint first and forward second makes
any leak visible: the forward run would report the strict estimator instead
of re-resolving to the adaptive one. This test runs that sequence through
openmc.lib in one process and checks the reported estimator each solve."""

import openmc
import openmc.lib
from openmc.examples import random_ray_three_region_cube


def test_random_ray_default_estimator_persistence(run_in_tmpdir, capfd):
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.particles = 10
    model.settings.inactive = 2
    model.settings.batches = 4
    # No volume_estimator set: both runs exercise the built-in default

    reported = []
    for adjoint in (True, False):
        model.settings.random_ray['adjoint'] = adjoint
        model.export_to_model_xml()
        openmc.lib.init()
        openmc.lib.run_random_ray()
        openmc.lib.finalize()
        out = capfd.readouterr().out
        for line in out.splitlines():
            if 'Volume Estimator Type' in line:
                reported.append(line.split('=')[-1].strip())

    # The adjoint run reports once per solve (forward-for-adjoint, adjoint);
    # the forward run reports once.
    assert reported == ['Strict Adaptive (auto)', 'Strict Adaptive (auto)',
                        'Adaptive (auto)'], (
        f"default volume estimator did not re-resolve across in-process "
        f"reruns: {reported}")
