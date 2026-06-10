import copy
import os

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


def configure_random_ray(model, mesh):
    """Convert the model to random ray mode, with settings shared by the
    weight window generation solve and the final solve."""
    model.convert_to_random_ray()
    rr = model.settings.random_ray
    rr['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]
    rr['distance_inactive'] = 100.0
    rr['distance_active'] = 400.0
    # convert_to_random_ray samples ray origins from the bounding box, which
    # for a spherical geometry includes regions outside the geometry. Sample
    # from the box inscribed in the sphere instead.
    ll = np.asarray(mesh.lower_left)
    ur = np.asarray(mesh.upper_right)
    center = (ll + ur) / 2
    half = (ur - ll).min() / 2 / np.sqrt(3.0) * 0.99
    rr['ray_source'] = openmc.IndependentSource(
        space=openmc.stats.Box(center - half, center + half))
    rr['source_shape'] = 'flat'
    rr['sample_method'] = 'halton'
    rr['volume_estimator'] = 'hybrid'
    model.settings.particles = 1000
    model.settings.batches = 20
    model.settings.inactive = 15


def test_random_ray_auto_convert_bootstrap():
    # Tests the weight window bootstrapped MGXS generation workflow: cheap
    # global weight windows (stochastic_slab MGXS + random ray FW-CADIS) feed a
    # material_wise generation via the weight_windows_file argument. The
    # example's steel pocket sits behind ~1 m of concrete, which the
    # material_wise generation solve essentially never reaches without weight
    # windows, so the pocket flux tallied in the final solve depends on the
    # bootstrapped steel cross sections being generated correctly.
    openmc.reset_auto_ids()
    model = sphere_with_shielded_pocket()

    # Inherited by the continuous energy generation solves below: weight window
    # checkpoints help the bootstrapped solve reach the pocket, and must have
    # no effect on the random ray solves (random ray generates weight windows
    # but never applies them).
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
    mg_model = copy.deepcopy(model)
    mg_model.convert_to_multigroup(
        method='stochastic_slab', groups=GROUPS, nparticles=50,
        overwrite_mgxs_library=True, mgxs_path='mgxs.h5')
    configure_random_ray(mg_model, mesh)
    mg_model.settings.weight_window_generators = openmc.WeightWindowGenerator(
        method='fw_cadis', mesh=mesh, max_realizations=20,
        energy_bounds=openmc.mgxs.EnergyGroups(GROUPS).group_edges)
    mg_model.run()

    # Bootstrap the material-wise MGXS generation with those weight windows
    model.convert_to_multigroup(
        method='material_wise', groups=GROUPS, nparticles=1,
        overwrite_mgxs_library=True, mgxs_path='mgxs.h5',
        weight_windows_file='weight_windows.h5')

    # Run a random ray solve with the bootstrapped library, tallying the flux
    # in every region
    configure_random_ray(model, mesh)
    tally = openmc.Tally(name='flux')
    tally.filters = [
        openmc.CellFilter(list(model.geometry.get_all_cells().values())),
        openmc.EnergyFilter(openmc.mgxs.EnergyGroups(GROUPS).group_edges)]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    harness = MGXSTestHarness('statepoint.20.h5', model)
    harness.main()
