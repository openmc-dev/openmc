import copy

import numpy as np
import openmc
from openmc.examples import sphere_with_shielded_pocket
from openmc.utility_funcs import change_directory


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


def make_weight_windows(model):
    """Generate FW-CADIS weight windows for the model with random ray, using
    MGXS from the stochastic_slab method, and return the weight window file."""
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
    mg_model.run()
    return 'weight_windows.h5'


def test_random_ray_auto_convert_bootstrap(tmp_path):
    """Bootstrapped material-wise MGXS generation with weight windows.

    The example model holds a steel pocket behind ~1 m of concrete, which a
    ~200 history analog solve reaches with probability well below 1%. The
    "material_wise" MGXS method therefore produces zero cross sections for the
    steel unless the continuous energy solve applies weight windows, supplied
    here via the weight_windows_file argument from a stochastic_slab + random
    ray FW-CADIS generation pass (the analog stage doubles as the negative
    control).
    """
    with change_directory(tmp_path):
        openmc.reset_auto_ids()
        model = sphere_with_shielded_pocket()

        # Inherited by every continuous energy solve below. Weight window
        # checkpoints assist the bootstrapped solve; they also propagate into
        # the random ray weight window generation run, where they must have no
        # effect on the rays (random ray generates weight windows but must
        # never apply them).
        model.settings.weight_window_checkpoints = {
            'collision': True, 'surface': True}
        model.settings.max_history_splits = 100_000

        # Stage 1: cheap global weight windows (stochastic_slab -> FW-CADIS)
        ww_file = make_weight_windows(model)

        # Stage 2 (negative control): the analog material_wise solve cannot
        # reach the steel pocket, so its MGXS must come out zero.
        analog = copy.deepcopy(model)
        analog.convert_to_multigroup(
            method='material_wise', groups=GROUPS, nparticles=1,
            overwrite_mgxs_library=True, mgxs_path='mgxs_analog.h5')
        covered, uncovered = mgxs_coverage('mgxs_analog.h5')
        assert uncovered == ['Steel']
        assert sorted(covered) == ['Air', 'Concrete']

        # Stage 3: the same solve bootstrapped with weight windows reaches the
        # pocket and produces a complete library.
        boot = copy.deepcopy(model)
        boot.convert_to_multigroup(
            method='material_wise', groups=GROUPS, nparticles=1,
            overwrite_mgxs_library=True, mgxs_path='mgxs_boot.h5',
            weight_windows_file=ww_file)
        covered, uncovered = mgxs_coverage('mgxs_boot.h5')
        assert uncovered == []
        assert sorted(covered) == ['Air', 'Concrete', 'Steel']
