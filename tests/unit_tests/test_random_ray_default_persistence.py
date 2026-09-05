"""The random ray volume-estimator setting must not leak across an
in-process finalize/re-initialize cycle (openmc.lib workflows such as
iterative weight window generation). The solver never modifies the
configured setting, as "auto" is resolved into a separate run-scoped value.
However, the configured setting is a static that XML parsing only assigns
when the element is present, so openmc_finalize_random_ray() must restore
the "auto" default between runs. Running an explicit estimator first and a
default model second makes a missed restore visible, as the second run
would report the first run's estimator instead of resolving "auto". The
default second run is an adjoint solve, which also pins the auto routing
under openmc.lib (it must resolve to the strict adaptive estimator, once
per solve)."""

import openmc
import openmc.lib
from openmc.examples import random_ray_three_region_cube


def test_random_ray_default_estimator_persistence(run_in_tmpdir, capfd):
    openmc.reset_auto_ids()
    model = random_ray_three_region_cube()
    model.settings.particles = 10
    model.settings.inactive = 2
    model.settings.batches = 4

    reported = []
    for explicit in (True, False):
        if explicit:
            model.settings.random_ray['volume_estimator'] = 'naive'
        else:
            del model.settings.random_ray['volume_estimator']
            model.settings.random_ray['adjoint'] = True
        model.export_to_model_xml()
        openmc.lib.init()
        openmc.lib.run_random_ray()
        openmc.lib.finalize()
        out = capfd.readouterr().out
        for line in out.splitlines():
            if 'Volume Estimator Type' in line:
                reported.append(line.split('=')[-1].strip())

    # The forward run reports once, while the adjoint run reports once per
    # solve (forward-for-adjoint, adjoint).
    assert reported == ['Naive', 'Strict Adaptive (auto)',
                        'Strict Adaptive (auto)'], (
        f"volume estimator leaked across in-process reruns: {reported}")
