import copy
import os

import numpy as np
import openmc
from openmc.examples import sphere_with_shielded_pocket

from tests.testing_harness import TolerantPyAPITestHarness


GROUPS = 'CASMO-4'


def mgxs_coverage(path):
    """Return (covered, uncovered) material name lists for an MGXS library.

    A material is "covered" if its total cross section is finite with at least
    one positive value; materials that received no tallies during generation
    have all-zero cross sections.
    """
    library = openmc.MGXSLibrary.from_hdf5(path)
    covered, uncovered = [], []
    for xsdata in library.xsdatas:
        total = np.asarray(xsdata.total[0], dtype=float)
        if np.all(np.isfinite(total)) and np.any(total > 0.0):
            covered.append(xsdata.name)
        else:
            uncovered.append(xsdata.name)
    return covered, uncovered


class MGXSBootstrapTestHarness(TolerantPyAPITestHarness):
    """Runs the full MGXS bootstrap workflow.

    The harness model is the random ray FW-CADIS weight window generation run
    (stage 1). After it executes, the analog (stage 2) and bootstrapped
    (stage 3) material-wise MGXS generation solves are run from the original
    continuous energy model.

    The compared results contain the generated weight window bounds and the
    analog MGXS data, both of which are reproducible across runs and thread
    counts, so they detect subtle changes anywhere in the random ray solver,
    the weight window generation machinery, or the MGXS generation chain. The
    bootstrapped MGXS values are deliberately NOT value-compared: a Monte
    Carlo solve with weight window splitting is not reproducible run-to-run
    when multithreaded (the splitting order perturbs the tracking), so only
    its binary material coverage -- which is the feature guarantee -- is
    asserted.
    """

    def __init__(self, statepoint_name, model, ce_model):
        super().__init__(statepoint_name, model)
        self._ce_model = ce_model
        self._analog_uncovered = None
        self._boot_uncovered = None

    def _run_openmc(self):
        # Stage 1: random ray FW-CADIS solve producing weight_windows.h5
        super()._run_openmc()

        # Stage 2 (negative control): the analog material_wise solve cannot
        # reach the steel pocket, so its MGXS must come out zero.
        analog = copy.deepcopy(self._ce_model)
        analog.convert_to_multigroup(
            method='material_wise', groups=GROUPS, nparticles=1,
            overwrite_mgxs_library=True, mgxs_path='mgxs_analog.h5')
        _, self._analog_uncovered = mgxs_coverage('mgxs_analog.h5')

        # Stage 3: the same solve bootstrapped with the stage 1 weight windows
        # reaches the pocket and produces a complete library.
        boot = copy.deepcopy(self._ce_model)
        boot.convert_to_multigroup(
            method='material_wise', groups=GROUPS, nparticles=1,
            overwrite_mgxs_library=True, mgxs_path='mgxs_boot.h5',
            weight_windows_file='weight_windows.h5')
        _, self._boot_uncovered = mgxs_coverage('mgxs_boot.h5')

        # Hard semantic asserts (rather than reference comparisons) so that a
        # broken feature can never be baked into the reference files by an
        # --update run.
        assert self._analog_uncovered == ['Steel']
        assert self._boot_uncovered == []

    def _get_results(self):
        """Digest the weight window bounds, the coverage verdicts, and the
        analog MGXS data."""
        # Weight window bounds from the random ray FW-CADIS generation
        ww = openmc.WeightWindowsList.from_hdf5()[0]
        lines = [str(ww.mesh), 'Lower Bounds']
        lines += [f'{x:.2e}' for x in ww.lower_ww_bounds.flatten()]
        lines.append('Upper Bounds')
        lines += [f'{x:.2e}' for x in ww.upper_ww_bounds.flatten()]

        # Material coverage verdicts from the two material_wise solves
        lines.append(f'analog uncovered: {sorted(self._analog_uncovered)}')
        lines.append(f'bootstrap uncovered: {sorted(self._boot_uncovered)}')

        # Analog MGXS data (the bootstrapped library is excluded; see class
        # docstring)
        library = openmc.MGXSLibrary.from_hdf5('mgxs_analog.h5')
        for xsdata in sorted(library.xsdatas, key=lambda x: x.name):
            lines.append(f'analog MGXS {xsdata.name}')
            for label, array in (('total', xsdata.total[0]),
                                 ('absorption', xsdata.absorption[0]),
                                 ('scatter', xsdata.scatter_matrix[0])):
                lines.append(label)
                lines += [f'{v:.6e}' for v in
                          np.asarray(array, dtype=float).flatten()]
        return '\n'.join(lines) + '\n'

    def _cleanup(self):
        super()._cleanup()
        for f in ('weight_windows.h5', 'mgxs_slab.h5', 'mgxs_analog.h5',
                  'mgxs_boot.h5'):
            if os.path.exists(f):
                os.remove(f)


def test_random_ray_auto_convert_bootstrap():
    """Bootstrapped material-wise MGXS generation with weight windows.

    The example model holds a steel pocket behind ~1 m of concrete, which a
    ~200 history analog solve reaches with probability well below 1%. The
    "material_wise" MGXS method therefore produces zero cross sections for the
    steel unless the continuous energy solve applies weight windows, supplied
    via the weight_windows_file argument from a stochastic_slab + random ray
    FW-CADIS generation pass.
    """
    openmc.reset_auto_ids()
    model = sphere_with_shielded_pocket()

    # Inherited by every continuous energy solve below. Weight window
    # checkpoints assist the bootstrapped solve; they also propagate into the
    # random ray weight window generation run, where they must have no effect
    # on the rays (random ray generates weight windows but never applies
    # them).
    model.settings.weight_window_checkpoints = {
        'collision': True, 'surface': True}
    model.settings.max_history_splits = 100_000

    # Convert a copy to multigroup with the stochastic_slab method and set up
    # the random ray FW-CADIS weight window generation run (stage 1).
    mg_model = copy.deepcopy(model)
    mg_model.convert_to_multigroup(
        method='stochastic_slab', groups=GROUPS, nparticles=50,
        overwrite_mgxs_library=True, mgxs_path='mgxs_slab.h5')
    mg_model.convert_to_random_ray()

    # Overlay a ~10 cm source region / weight window mesh on the geometry
    bbox = mg_model.geometry.bounding_box
    ll = np.asarray(bbox.lower_left)
    ur = np.asarray(bbox.upper_right)
    mesh = openmc.RegularMesh()
    mesh.dimension = np.round((ur - ll) / 10.0).astype(int)
    mesh.lower_left = ll
    mesh.upper_right = ur

    rr = mg_model.settings.random_ray
    rr['source_region_meshes'] = [(mesh, [mg_model.geometry.root_universe])]
    rr['distance_inactive'] = 100.0
    rr['distance_active'] = 400.0
    # convert_to_random_ray samples ray origins from the bounding box, which
    # for a spherical geometry includes regions outside the geometry. Sample
    # from the box inscribed in the sphere instead.
    center = (ll + ur) / 2
    half = (ur - ll).min() / 2 / np.sqrt(3.0) * 0.99
    rr['ray_source'] = openmc.IndependentSource(
        space=openmc.stats.Box(center - half, center + half))
    rr['source_shape'] = 'flat'
    rr['sample_method'] = 'halton'
    rr['volume_estimator'] = 'hybrid'

    mg_model.settings.particles = 1000
    mg_model.settings.batches = 20
    mg_model.settings.inactive = 15
    mg_model.settings.weight_window_generators = openmc.WeightWindowGenerator(
        method='fw_cadis', mesh=mesh, max_realizations=20,
        energy_bounds=openmc.mgxs.EnergyGroups(GROUPS).group_edges)

    harness = MGXSBootstrapTestHarness('statepoint.20.h5', mg_model, model)
    harness.main()
