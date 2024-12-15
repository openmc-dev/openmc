import os
from pathlib import Path

import numpy as np
import pytest
from uncertainties import ufloat

import openmc
import openmc.lib
from openmc.stats import Discrete, Point

from tests import cdtemp


@pytest.fixture
def wws():

    # weight windows
    ww_files = ("ww_n.txt", "ww_p.txt")
    cwd = Path(__file__).parent.absolute()
    ww_n_file, ww_p_file = [cwd / Path(f) for f in ww_files]

    # load pre-generated weight windows
    # (created using the same tally as above)
    ww_n_lower_bnds = np.loadtxt(ww_n_file)
    ww_p_lower_bnds = np.loadtxt(ww_p_file)

    # create a mesh matching the one used
    # to generate the weight windows
    ww_mesh = openmc.RegularMesh()
    ww_mesh.lower_left = (-240, -240, -240)
    ww_mesh.upper_right = (240, 240, 240)
    ww_mesh.dimension = (5, 6, 7)

    # energy bounds matching those of the
    # generated weight windows
    e_bnds = [0.0, 0.5, 2e7]

    ww_n = openmc.WeightWindows(
        ww_mesh, ww_n_lower_bnds, None, 10.0, e_bnds, survival_ratio=1.01
    )

    ww_p = openmc.WeightWindows(
        ww_mesh, ww_p_lower_bnds, None, 10.0, e_bnds, survival_ratio=1.01
    )

    return [ww_n, ww_p]


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    # materials (M4 steel alloy)
    m4 = openmc.Material()
    m4.set_density("g/cc", 2.3)
    m4.add_nuclide("H1", 0.168018676)
    m4.add_nuclide("H2", 1.93244e-05)
    m4.add_nuclide("O16", 0.561814465)
    m4.add_nuclide("O17", 0.00021401)
    m4.add_nuclide("Na23", 0.021365)
    m4.add_nuclide("Al27", 0.021343)
    m4.add_nuclide("Si28", 0.187439342)
    m4.add_nuclide("Si29", 0.009517714)
    m4.add_nuclide("Si30", 0.006273944)
    m4.add_nuclide("Ca40", 0.018026179)
    m4.add_nuclide("Ca42", 0.00012031)
    m4.add_nuclide("Ca43", 2.51033e-05)
    m4.add_nuclide("Ca44", 0.000387892)
    m4.add_nuclide("Ca46", 7.438e-07)
    m4.add_nuclide("Ca48", 3.47727e-05)
    m4.add_nuclide("Fe54", 0.000248179)
    m4.add_nuclide("Fe56", 0.003895875)
    m4.add_nuclide("Fe57", 8.99727e-05)
    m4.add_nuclide("Fe58", 1.19737e-05)

    s0 = openmc.Sphere(r=240)
    s1 = openmc.Sphere(r=250, boundary_type="vacuum")

    c0 = openmc.Cell(fill=m4, region=-s0)
    c1 = openmc.Cell(region=+s0 & -s1)

    model.geometry = openmc.Geometry([c0, c1])

    # settings
    settings = model.settings
    settings.run_mode = "fixed source"
    settings.particles = 500
    settings.batches = 2
    settings.max_history_splits = 100
    settings.photon_transport = True
    space = Point((0.001, 0.001, 0.001))
    energy = Discrete([14e6], [1.0])

    settings.source = openmc.IndependentSource(space=space, energy=energy)

    # tally
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-240, -240, -240)
    mesh.upper_right = (240, 240, 240)
    mesh.dimension = (3, 5, 7)

    mesh_filter = openmc.MeshFilter(mesh)

    e_bnds = [0.0, 0.5, 2e7]
    energy_filter = openmc.EnergyFilter(e_bnds)

    particle_filter = openmc.ParticleFilter(["neutron", "photon"])

    tally = openmc.Tally()
    tally.filters = [mesh_filter, energy_filter, particle_filter]
    tally.scores = ["flux"]

    model.tallies.append(tally)

    return model


def test_weightwindows(model, wws):

    ww_files = ("ww_n.txt", "ww_p.txt")
    cwd = Path(__file__).parent.absolute()
    filepaths = [cwd / Path(f) for f in ww_files]

    with cdtemp(filepaths):
        # run once with variance reduction off
        model.settings.weight_windows_on = False
        analog_sp = model.run()
        os.rename(analog_sp, "statepoint.analog.h5")

        model.settings.weight_windows = wws

        # check that string form of the class can be created
        for ww in model.settings.weight_windows:
            str(ww)

        # run again with variance reduction on
        model.settings.weight_windows_on = True
        ww_sp = model.run()
        os.rename(ww_sp, "statepoint.ww.h5")

        # load both statepoints and examine results
        asp = openmc.StatePoint("statepoint.analog.h5")
        wsp = openmc.StatePoint("statepoint.ww.h5")

        analog_tally = asp.tallies[1]
        ww_tally = wsp.tallies[1]

        def compare_results(particle, analog_tally, ww_tally):
            # get values from each of the tallies
            an_mean = analog_tally.get_values(
                filters=[openmc.ParticleFilter], filter_bins=[(particle,)]
            )
            ww_mean = ww_tally.get_values(
                filters=[openmc.ParticleFilter], filter_bins=[(particle,)]
            )

            # expect that more bins were scored with weight windows than
            # the analog run
            assert np.count_nonzero(an_mean) < np.count_nonzero(ww_mean)

            an_rel_err = analog_tally.get_values(
                filters=[openmc.ParticleFilter],
                filter_bins=[(particle,)],
                value="rel_err",
            )
            ww_rel_err = ww_tally.get_values(
                filters=[openmc.ParticleFilter],
                filter_bins=[(particle,)],
                value="rel_err",
            )

            an_rel_err[an_mean == 0.0] = 1.0
            ww_rel_err[ww_mean == 0.0] = 1.0

            an_avg_rel_err = np.mean(an_rel_err)
            ww_avg_rel_err = np.mean(ww_rel_err)

            # expect that the average relative error in the tally
            # decreases
            assert an_avg_rel_err > ww_avg_rel_err

            # ensure that the value of the mesh bin containing the
            # source is statistically similar in both runs
            an_std_dev = analog_tally.get_values(
                filters=[openmc.ParticleFilter],
                filter_bins=[(particle,)],
                value="std_dev",
            )
            ww_std_dev = ww_tally.get_values(
                filters=[openmc.ParticleFilter],
                filter_bins=[(particle,)],
                value="std_dev",
            )

            # index of the mesh bin containing the source for the higher
            # energy group
            source_bin_idx = (an_mean.shape[0] // 2, 0, 0)

            an_source_bin = ufloat(an_mean[source_bin_idx], an_std_dev[source_bin_idx])
            ww_source_bin = ufloat(ww_mean[source_bin_idx], ww_std_dev[source_bin_idx])

            diff = an_source_bin - ww_source_bin

            # check that values are within two combined standard deviations
            assert abs(diff.nominal_value) / diff.std_dev < 2.0

        compare_results("neutron", analog_tally, ww_tally)
        compare_results("photon", analog_tally, ww_tally)


def test_lower_ww_bounds_shape():
    """checks that lower_ww_bounds is reshaped to the mesh dimension when set"""
    ww_mesh = openmc.RegularMesh()
    ww_mesh.lower_left = (-10, -10, -10)
    ww_mesh.upper_right = (10, 10, 10)
    ww_mesh.dimension = (2, 3, 4)

    ww = openmc.WeightWindows(
        mesh=ww_mesh,
        lower_ww_bounds=[1] * 24,
        upper_bound_ratio=5,
        energy_bounds=(1, 1e40),
    )
    assert ww.lower_ww_bounds.shape == (2, 3, 4, 1)


def test_roundtrip(run_in_tmpdir, model, wws):
    model.settings.weight_windows = wws

    # write the model with weight windows to XML
    model.export_to_xml()

    # ensure that they can be read successfully from XML and that they match the input values
    model_read = openmc.Model.from_xml()

    zipped_wws = zip(model.settings.weight_windows, model_read.settings.weight_windows)

    # ensure the lower bounds read in from the XML match those of the
    for ww_out, ww_in in zipped_wws:
        assert ww_out == ww_in


def test_ww_attrs_python(model):
    mesh = openmc.RegularMesh.from_domain(model.geometry)
    lower_bounds = np.ones(mesh.dimension)

    # ensure that creation of weight window objects with default arg values
    # is successful
    wws = openmc.WeightWindows(mesh, lower_bounds, upper_bound_ratio=10.0)

    assert wws.energy_bounds == None

    wwg = openmc.WeightWindowGenerator(mesh)

    assert wwg.energy_bounds == None


def test_ww_attrs_capi(run_in_tmpdir, model):
    model.export_to_xml()

    openmc.lib.init()

    tally = openmc.lib.tallies[model.tallies[0].id]

    wws = openmc.lib.WeightWindows.from_tally(tally)

    # this is the first weight window object created
    assert wws.id == 1

    with pytest.raises(ValueError):
        tally.find_filter(openmc.lib.AzimuthalFilter)

    mesh_filter = tally.find_filter(openmc.lib.MeshFilter)
    mesh = mesh_filter.mesh

    assert wws.mesh.id == mesh.id

    assert wws.particle == openmc.ParticleType.NEUTRON

    wws.particle = 1
    assert wws.particle == openmc.ParticleType.PHOTON
    wws.particle = "photon"
    assert wws.particle == openmc.ParticleType.PHOTON

    with pytest.raises(ValueError):
        wws.particle = "🌠"

    energy_filter = tally.find_filter(openmc.lib.EnergyFilter)
    np.testing.assert_allclose(np.unique(energy_filter.bins), wws.energy_bounds)

    # at this point the weight window bounds are uninitialized
    assert all(wws.bounds[0] == -1)
    assert all(wws.bounds[1] == -1)

    wws = openmc.lib.WeightWindows.from_tally(tally, particle="photon")
    assert wws.id == 2
    assert wws.particle == openmc.ParticleType.PHOTON

    openmc.lib.finalize()


@pytest.mark.parametrize("library", ("libmesh", "moab"))
def test_unstructured_mesh_applied_wws(request, run_in_tmpdir, library):
    """
    Ensure that weight windows on unstructured mesh work when
    they aren't part of a tally or weight window generator
    """

    if library == "libmesh" and not openmc.lib._libmesh_enabled():
        pytest.skip("LibMesh not enabled in this build.")
    if library == "moab" and not openmc.lib._dagmc_enabled():
        pytest.skip("DAGMC (and MOAB) mesh not enabled in this build.")

    water = openmc.Material(name="water")
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cc", 1.0)
    box = openmc.model.RectangularParallelepiped(
        *(3 * [-10, 10]), boundary_type="vacuum"
    )
    cell = openmc.Cell(region=-box, fill=water)

    geometry = openmc.Geometry([cell])
    mesh_file = str(request.fspath.dirpath() / "test_mesh_tets.exo")
    mesh = openmc.UnstructuredMesh(mesh_file, library)

    dummy_wws = np.ones((12_000,))

    wws = openmc.WeightWindows(mesh, dummy_wws, upper_bound_ratio=5.0)

    model = openmc.Model(geometry)
    model.settings.weight_windows = wws
    model.settings.weight_windows_on = True
    model.settings.run_mode = "fixed source"
    model.settings.particles = 100
    model.settings.batches = 2
    model.run()
