import copy
import os
from pathlib import Path

import numpy as np
import openmc
import openmc.mgxs
from openmc.examples import sphere_with_shielded_pocket

from tests.testing_harness import TolerantPyAPITestHarness


GROUPS = 'CASMO-4'


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        for f in ('mgxs.h5', 'weight_windows.h5'):
            if os.path.exists(f):
                os.remove(f)


def generate_weight_windows(model, mesh):
    """Convert a multigroup model to random ray and run a FW-CADIS weight
    window generation solve, producing weight_windows.h5."""
    model.convert_to_random_ray()
    rr = model.settings.random_ray
    rr['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]
    rr['distance_inactive'] = 100.0
    rr['distance_active'] = 400.0
    rr['source_shape'] = 'flat'
    rr['sample_method'] = 'halton'
    rr['volume_estimator'] = 'hybrid'
    model.settings.particles = 1000
    model.settings.batches = 20
    model.settings.inactive = 15
    model.settings.weight_window_generators = openmc.WeightWindowGenerator(
        method='fw_cadis', mesh=mesh, max_realizations=20,
        energy_bounds=openmc.mgxs.EnergyGroups(GROUPS).group_edges)
    model.run()


def test_random_ray_auto_convert_bootstrap():
    # Tests the weight window bootstrapped MGXS generation workflow, which
    # consists of five OpenMC runs: a stochastic_slab MGXS generation, a random
    # ray FW-CADIS weight window generation, a material_wise MGXS generation
    # using those weight windows (which reach the steel pocket that an analog
    # solve essentially never tallies), a second FW-CADIS generation using the
    # bootstrapped library, and a final continuous energy Monte Carlo solve
    # using the improved weight windows, whose tallies are compared against
    # the reference results.
    openmc.reset_auto_ids()
    model = sphere_with_shielded_pocket()

    # Used by the continuous energy solves below; these have no effect on the
    # random ray solves (random ray generates weight windows but never applies
    # them).
    model.settings.weight_window_checkpoints = {
        'collision': True, 'surface': True}
    model.settings.max_history_splits = 100_000

    # Overlay a ~10 cm source region / weight window mesh on the geometry
    bbox = model.geometry.bounding_box
    mesh = openmc.RegularMesh()
    mesh.dimension = np.round(
        (np.asarray(bbox.upper_right) - np.asarray(bbox.lower_left)) / 10.0
    ).astype(int)
    mesh.lower_left = bbox.lower_left
    mesh.upper_right = bbox.upper_right

    # Generate weight windows covering the whole problem with the
    # stochastic_slab method and a random ray FW-CADIS solve
    slab_model = copy.deepcopy(model)
    slab_model.convert_to_multigroup(
        method='stochastic_slab', groups=GROUPS, particles=50,
        overwrite_mgxs_library=True, mgxs_path='mgxs.h5')
    generate_weight_windows(slab_model, mesh)

    # Bootstrap the material-wise MGXS generation with those weight windows,
    # then regenerate the weight windows from the higher-fidelity library
    boot_model = copy.deepcopy(model)
    boot_model.convert_to_multigroup(
        groups=GROUPS, particles=1,
        weight_windows_file=Path('weight_windows.h5').resolve(),
        overwrite_mgxs_library=True, mgxs_path='mgxs.h5')
    generate_weight_windows(boot_model, mesh)

    # Run the continuous energy model with the improved weight windows,
    # tallying the flux in every region
    model.settings.weight_windows_file = 'weight_windows.h5'
    model.settings.particles = 20
    model.settings.batches = 10
    tally = openmc.Tally(name='flux')
    tally.filters = [
        openmc.CellFilter(list(model.geometry.get_all_cells().values()))]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
