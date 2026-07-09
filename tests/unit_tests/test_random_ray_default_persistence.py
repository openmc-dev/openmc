"""The random ray volume-estimator default must survive an in-process
finalize/re-initialize cycle (openmc.lib workflows such as iterative weight
window generation). openmc_finalize_random_ray() restores the built-in
defaults between runs; if it restores a different estimator than the static
default, the first and subsequent runs of a process silently use different
estimators. This test runs the same model twice through openmc.lib in one
process and checks the reported estimator both times."""

import openmc
import openmc.lib
from openmc.examples import random_ray_three_region_cube


def test_random_ray_default_estimator_persistence(run_in_tmpdir, capfd):
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.particles = 10
    model.settings.inactive = 2
    model.settings.batches = 4
    # No volume_estimator set: exercises the built-in default both runs
    model.export_to_model_xml()

    reported = []
    for _ in range(2):
        openmc.lib.init()
        openmc.lib.run_random_ray()
        openmc.lib.finalize()
        out = capfd.readouterr().out
        for line in out.splitlines():
            if 'Volume Estimator Type' in line:
                reported.append(line.split('=')[-1].strip())

    assert reported == ['Adaptive', 'Adaptive'], (
        f"default volume estimator changed across in-process reruns: "
        f"{reported}")
