import copy
from itertools import permutations
from pathlib import Path

import numpy as np
import openmc
import openmc.lib
import pytest


@pytest.fixture
def model():
    openmc.reset_auto_ids()

    # create a simple spherical shell shielding model

    ### Materials ###
    water = openmc.Material()
    water.set_density('g/cc', 1.0)
    water.add_nuclide("H1", 2)
    water.add_nuclide("O16", 1)

    steel = openmc.Material()
    steel.set_density('g/cc', 8.0)
    steel.add_nuclide("Fe56", 1.0)

    air = openmc.Material()
    air.set_density('g/cc', 0.001205)
    air.add_nuclide("N14", 0.781557629247)
    air.add_nuclide("O16", 0.210668126508)

    boron = openmc.Material()
    boron.set_density('g/cc', 2.52)
    boron.add_nuclide("B10", 0.15856)
    boron.add_nuclide("B11", 0.64144)
    boron.add_nuclide("C0", 0.2)

    ### Geometry ###
    radii = [5.0, 10.0, 30.0, 31.0, 50.0]

    surfs = [openmc.Sphere(r=r) for r in radii]

    surfs[-1].boundary_type = 'vacuum'

    regions = openmc.model.subdivide(surfs)

    mats = [air, water, steel, boron, air]

    cells = [openmc.Cell(fill=m, region=r) for r, m in zip(regions, mats)]

    geometry = openmc.Geometry(cells)

    ### Settings ###

    settings = openmc.Settings(
        run_mode='fixed source',
        particles=100,
        batches=10,
        max_history_splits=10,
        survival_biasing=False
    )

    # 10 keV neutron point source at the origin
    space = openmc.stats.Point()
    energy = openmc.stats.Discrete(x=[1e4], p=[1.0])
    settings.source = openmc.IndependentSource(space=space, energy=energy)

    return openmc.Model(geometry=geometry, settings=settings)


# create a tally used for weight window generation
mesh = openmc.RegularMesh()
mesh.lower_left = [-50.0] * 3
mesh.upper_right = [50.0] * 3
# use un-equal mesh widths in each dimension to more robustly check
# use of tally data
mesh.dimension = (19, 20, 21)

mf = openmc.MeshFilter(mesh)

ef = openmc.EnergyFilter([0.0, 1e7])

pf = openmc.ParticleFilter(['neutron', 'photon'])

filters = [mf, ef, pf]

test_cases = list(permutations(filters))
test_cases += list(permutations(filters[:-1]))
test_cases += list(permutations(filters[::2]))


def labels(params):
    out = []
    for p in params:
        if isinstance(p, openmc.ParticleFilter):
            out.append('particle')
        elif isinstance(p, openmc.MeshFilter):
            out.append('mesh')
        elif isinstance(p, openmc.EnergyFilter):
            out.append('energy')
    return "filters:" + '-'.join(out)


@pytest.mark.parametrize("filters", test_cases, ids=labels)
def test_ww_gen(filters, run_in_tmpdir, model):

    tally = openmc.Tally()
    tally.filters = list(filters)
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    model.export_to_model_xml()

    ref_lower = None
    ref_upper = None
    # test weight window generation capability
    with openmc.lib.run_in_memory():

        # retrieve the tally we created above in memory
        lib_tally = openmc.lib.tallies[tally.id]

        # create a new weight windows object
        ww = openmc.lib.WeightWindows.from_tally(lib_tally)

        # run particle transport
        openmc.lib.run()

        # capture analog data
        analog_mean = np.copy(lib_tally.mean)

        # update the weight window values using tally results
        ww.update_magic(lib_tally)

        assert any(ww.bounds[0] != -1)
        assert any(ww.bounds[1] != -1)

        # make sure that the weight window update doesn't change tally values
        np.testing.assert_equal(lib_tally.mean, analog_mean)

        # check against weight windows from the previous iteration
        # the order of filters should not change the weight window values
        if ref_lower is None:
            ref_lower = ww.bounds[0].copy()
        else:
            np.testing.assert_equal(ref_lower, ww.bounds[0])

        if ref_upper is None:
            ref_upper = ww.bounds[1].copy()
        else:
            np.testing.assert_equal(ref_upper, ww.bounds[1])

        # turn on weight windows for the subsequent run
        openmc.lib.settings.weight_windows_on = True

        openmc.lib.hard_reset()

        openmc.lib.run()

        ww_mean = np.copy(lib_tally.mean)

        # we expect that the application of weight windows will populate more tally
        # bins than the analog run for the meshes in this test model
        assert any(ww_mean != analog_mean)
        assert np.count_nonzero(ww_mean) > np.count_nonzero(analog_mean)


def test_ww_import_export(run_in_tmpdir, model):
    # create a tally for weight windows
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-50.0] * 3
    mesh.upper_right = [50.0] * 3
    # use un-equal mesh widths in each dimension to more robustly check
    # use of tally data
    mesh.dimension = (3, 4, 5)

    mf = openmc.MeshFilter(mesh)

    e_groups = np.logspace(0, 7, 8)
    ef = openmc.EnergyFilter(e_groups)

    pf = openmc.ParticleFilter(['neutron', 'photon'])

    tally = openmc.Tally()
    tally.filters = [mf, ef, pf]
    tally.scores = ['flux']

    model.tallies = openmc.Tallies([tally])

    # first, generate some weight windows
    model.export_to_model_xml()

    openmc.lib.init()

    tally = openmc.lib.tallies[tally.id]

    ww = openmc.lib.WeightWindows.from_tally(tally)

    openmc.lib.run()

    mean_before = np.array(tally.mean)

    ww.update_magic(tally)

    mean_after = np.array(tally.mean)

    assert (mean_before == mean_after).all()

    lb_before, up_before = ww.bounds

    # set some additional weight windows properties after transport
    ww.survival_ratio = 0.7
    assert ww.survival_ratio == 0.7

    ww.weight_cutoff = 1e-10
    assert ww.weight_cutoff == 1e-10

    ww.max_lower_bound_ratio = 200.0
    assert ww.max_lower_bound_ratio == 200.0

    ww.max_split = 26000
    assert ww.max_split == 26000

    openmc.lib.export_weight_windows()

    assert Path('weight_windows.h5').exists()

    openmc.lib.import_weight_windows('weight_windows.h5')

    imported_ww = openmc.lib.weight_windows[2]

    lb_after, up_after = imported_ww.bounds

    assert np.allclose(lb_before, lb_after)
    assert np.allclose(up_before, up_after)

    assert ww.survival_ratio == imported_ww.survival_ratio
    assert ww.max_lower_bound_ratio == imported_ww.max_lower_bound_ratio
    assert ww.weight_cutoff == imported_ww.weight_cutoff
    assert ww.max_split == imported_ww.max_split

    openmc.lib.finalize()


def test_ww_gen_roundtrip(run_in_tmpdir, model):

    mesh = openmc.RegularMesh.from_domain(model.geometry.root_universe)
    energy_bounds = np.linspace(0.0, 1e8, 11)
    particle_type = 'neutron'

    wwg = openmc.WeightWindowGenerator(mesh, energy_bounds, particle_type)
    wwg.update_parameters = {'ratio' : 5.0,
                             'threshold': 0.8,
                             'value' : 'mean'}

    model.settings.weight_window_generators = wwg
    model.export_to_xml()

    model_in = openmc.Model.from_xml()

    assert len(model_in.settings.weight_window_generators) == 1

    wwg_in = model.settings.weight_window_generators[0]

    # rountrip tests
    model_in = openmc.Model.from_xml()
    assert len(model_in.settings.weight_window_generators) == 1
    wwg_in = model_in.settings.weight_window_generators[0]
    assert wwg_in.max_realizations == 1
    assert wwg_in.on_the_fly == True
    assert wwg_in.update_interval == 1
    assert wwg_in.update_parameters == wwg.update_parameters

    with pytest.raises(ValueError):
        wwg.method = '🦍🐒'

    with pytest.raises(TypeError):
        wwg.update_parameters = {'ratio' : 'one-to-one'}

    with pytest.raises(ValueError):
        wwg.max_realizations = -1


def test_python_hdf5_roundtrip(run_in_tmpdir, model):

    # add a tally to the model
    mesh = openmc.RegularMesh.from_domain(model.geometry)

    # some arbitrary energy groups
    e_groups = np.logspace(0, 6, 4)
    energy_filter = openmc.EnergyFilter(e_groups)

    bounds = np.arange(energy_filter.num_bins * np.prod(mesh.dimension))

    wws = openmc.WeightWindows(mesh, bounds, bounds, energy_bounds=e_groups)

    model.settings.weight_windows = [wws]

    model.export_to_xml()

    # initialize and export wws to HDF5
    openmc.lib.init()

    openmc.lib.export_weight_windows()

    openmc.lib.finalize()

    wws_hdf5 = openmc.WeightWindowsList.from_hdf5()[0]

    # ensure
    assert all(wws.energy_bounds == wws_hdf5.energy_bounds)
    assert wws.id == wws_hdf5.id
    assert wws.mesh.id == wws_hdf5.mesh.id
    np.testing.assert_array_equal(wws.lower_ww_bounds, wws_hdf5.lower_ww_bounds)
    np.testing.assert_array_equal(wws.upper_ww_bounds, wws_hdf5.upper_ww_bounds)


def test_ww_bounds_set_in_memory(run_in_tmpdir, model):
    tally = openmc.Tally()
    tally.filters = filters
    tally.scores = ['flux']
    model.tallies = [tally]

    bounds = np.arange(ef.num_bins * np.prod(mf.mesh.dimension))

    model.export_to_xml()

    openmc.lib.init()

    lib_tally = openmc.lib.tallies[tally.id]

    wws = openmc.lib.WeightWindows.from_tally(lib_tally)

    wws.bounds = (bounds, bounds)

    openmc.lib.finalize()


@pytest.mark.skipif(not openmc.lib._dagmc_enabled(), reason="DAGMC CAD geometry is not enabled.")
def test_ww_generation_with_dagmc(run_in_tmpdir):
    mat1 = openmc.Material(name="1")
    mat1.add_nuclide("H1", 1, percent_type="ao")
    mat1.set_density("g/cm3", 0.001)

    materials = openmc.Materials([mat1])
    dag_univ = openmc.DAGMCUniverse(
        Path(__file__).parent.parent / "dagmc" / "dagmc_tetrahedral_no_graveyard.h5m")
    bound_dag_univ = dag_univ.bounded_universe(padding_distance=1)
    geometry = openmc.Geometry(bound_dag_univ)

    settings = openmc.Settings()
    settings.batches = 6
    settings.particles = 30
    settings.run_mode = "fixed source"

    # Create a point source which are supported by random ray mode
    my_source = openmc.IndependentSource()
    my_source.space = openmc.stats.Point((0.25, 0.25, 0.25))
    my_source.energy = openmc.stats.delta_function(14e6)
    settings.source = my_source

    model = openmc.Model(geometry, materials, settings)

    rr_model = copy.deepcopy(model)
    rr_model.settings.inactive = 3

    rr_model.convert_to_multigroup(
        method="stochastic_slab",
        overwrite_mgxs_library=True,
        nparticles=10,
        groups="CASMO-2"
    )

    rr_model.convert_to_random_ray()

    mesh = openmc.RegularMesh.from_domain(rr_model, dimension=(4, 4, 4))

    # avoid writing files we don't make use of
    rr_model.settings.output = {"summary": False, "tallies": False}

    # Subdivide random ray source regions
    rr_model.settings.random_ray["source_region_meshes"] = [
        (mesh, [rr_model.geometry.root_universe])
    ]

    # less likely to get negative values in the weight window
    rr_model.settings.random_ray["volume_estimator"] = "naive"

    # Add a weight window generator to the model
    rr_model.settings.weight_window_generators = openmc.WeightWindowGenerator(
        method="fw_cadis",
        mesh=mesh,
        max_realizations=42,
        particle_type='neutron',
        energy_bounds=[0.0, 100e6]
    )

    rr_model.run()

    model.settings.weight_windows_on = True
    model.settings.weight_window_checkpoints = {"collision": True, "surface": True}
    model.settings.survival_biasing = False
    model.settings.weight_windows = openmc.WeightWindowsList.from_hdf5()

    model.run()
