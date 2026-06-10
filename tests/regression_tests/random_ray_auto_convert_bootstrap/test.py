import copy
import os
from contextlib import contextmanager

import numpy as np
import openmc
from openmc.examples import sphere_with_shielded_pocket

from tests.regression_tests import config
from tests.testing_harness import TolerantPyAPITestHarness


GROUPS = 'CASMO-4'


@contextmanager
def single_thread():
    """Pin subprocess OpenMC runs to a single thread.

    The stages that produce the weight windows must be bit-reproducible:
    multithreaded tally and random ray flux accumulations carry last-ulp
    ordering jitter, and the weight window split/roulette decisions downstream
    amplify even one-ulp window perturbations into order-unity tally changes
    (particle weights are set from the window bounds, so weight-window
    comparisons can sit exactly on decision boundaries).
    """
    old = os.environ.get('OMP_NUM_THREADS')
    os.environ['OMP_NUM_THREADS'] = '1'
    try:
        yield
    finally:
        if old is None:
            del os.environ['OMP_NUM_THREADS']
        else:
            os.environ['OMP_NUM_THREADS'] = old


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

    The compared results contain the stochastic slab, analog, and bootstrapped
    MGXS libraries along with the generated weight window bounds, so subtle
    changes anywhere in the chain -- the MGXS generation methods, the random
    ray solver, the weight window generation, or the weight window application
    (splitting/roulette) in the Monte Carlo solve -- are detected. The whole
    chain is reproducible across runs and thread counts because the stages
    that produce the weight windows are pinned to a single thread (see
    single_thread); given bit-identical windows, the weight-window-split Monte
    Carlo solve is itself reproducible to the last ulp across thread counts.
    """

    def __init__(self, statepoint_name, model, ce_model):
        super().__init__(statepoint_name, model)
        self._ce_model = ce_model
        self._analog_uncovered = None
        self._boot_uncovered = None

    def _run_openmc(self):
        # Stage 1: random ray FW-CADIS solve producing weight_windows.h5,
        # single-threaded so the window bounds are bit-reproducible.
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(threads=1, openmc_exec=config['exe'],
                       mpi_args=mpi_args, event_based=config['event'])
        else:
            openmc.run(threads=1, openmc_exec=config['exe'],
                       event_based=config['event'])

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

        # Hard semantic asserts (in addition to the reference comparison) so
        # that a broken feature can never be baked into the reference files by
        # an --update run.
        assert self._analog_uncovered == ['Steel']
        assert self._boot_uncovered == []

    def _mgxs_lines(self, label, path):
        lines = []
        library = openmc.MGXSLibrary.from_hdf5(path)
        for xsdata in sorted(library.xsdatas, key=lambda x: x.name):
            lines.append(f'{label} MGXS {xsdata.name}')
            for name, array in (('total', xsdata.total[0]),
                                ('absorption', xsdata.absorption[0]),
                                ('scatter', xsdata.scatter_matrix[0])):
                lines.append(name)
                lines += [f'{v:.6e}' for v in
                          np.asarray(array, dtype=float).flatten()]
        return lines

    def _get_results(self):
        """Digest the coverage verdicts, the three MGXS libraries, and the
        weight window bounds."""
        lines = [f'analog uncovered: {sorted(self._analog_uncovered)}',
                 f'bootstrap uncovered: {sorted(self._boot_uncovered)}']
        lines += self._mgxs_lines('slab', 'mgxs_slab.h5')
        lines += self._mgxs_lines('analog', 'mgxs_analog.h5')
        lines += self._mgxs_lines('bootstrap', 'mgxs_boot.h5')

        ww = openmc.WeightWindowsList.from_hdf5()[0]
        lines.append(str(ww.mesh))
        lines.append('Lower Bounds')
        lines += [f'{x:.6e}' for x in ww.lower_ww_bounds.flatten()]
        lines.append('Upper Bounds')
        lines += [f'{x:.6e}' for x in ww.upper_ww_bounds.flatten()]
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
    # the random ray FW-CADIS weight window generation run (stage 1). The slab
    # solve feeds the weight windows, so it runs single-threaded to be
    # bit-reproducible (see single_thread).
    mg_model = copy.deepcopy(model)
    with single_thread():
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
