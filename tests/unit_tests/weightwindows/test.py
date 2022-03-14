import os
from pathlib import Path

import numpy as np
import pytest
from uncertainties import ufloat

import openmc
from openmc.stats import Discrete, Point

from tests import cdtemp

@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    # materials (M4 steel alloy)
    m4 = openmc.Material()
    m4.set_density('g/cc', 2.3)
    m4.add_nuclide('H1', 0.168018676)
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
    s1 = openmc.Sphere(r=250, boundary_type='vacuum')

    c0 = openmc.Cell(fill=m4, region=-s0)
    c1 = openmc.Cell(region=+s0 & -s1)

    model.geometry = openmc.Geometry([c0, c1])

    # settings
    settings = model.settings
    settings.run_mode = 'fixed source'
    settings.particles = 500
    settings.batches = 2
    settings.max_splits = 100
    settings.photon_transport = True
    space = Point((0.001, 0.001, 0.001))
    energy = Discrete([14E6], [1.0])

    settings.source = openmc.Source(space=space, energy=energy)

    # tally
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-240, -240, -240)
    mesh.upper_right = (240, 240, 240)
    mesh.dimension = (3, 5, 7)

    mesh_filter = openmc.MeshFilter(mesh)

    e_bnds = [0.0, 0.5, 2E7]
    energy_filter = openmc.EnergyFilter(e_bnds)

    particle_filter = openmc.ParticleFilter(['neutron', 'photon'])

    tally = openmc.Tally()
    tally.filters = [mesh_filter, energy_filter, particle_filter]
    tally.scores = ['flux']

    model.tallies.append(tally)

    return model


def test_weightwindows(model):

    ww_files = ('ww_n.txt', 'ww_p.txt')
    cwd = Path(__file__).parent.absolute()
    filepaths = [cwd / Path(f) for f in ww_files]

    with cdtemp(filepaths):
        # run once with variance reduction off
        model.settings.weight_windows_on = False
        analog_sp = model.run()
        os.rename(analog_sp, 'statepoint.analog.h5')

        # weight windows

        # load pre-generated weight windows
        # (created using the same tally as above)
        ww_n_lower_bnds = np.loadtxt('ww_n.txt')
        ww_p_lower_bnds = np.loadtxt('ww_p.txt')

        # create a mesh matching the one used
        # to generate the weight windows
        ww_mesh = openmc.RegularMesh()
        ww_mesh.lower_left = (-240, -240, -240)
        ww_mesh.upper_right = (240, 240, 240)
        ww_mesh.dimension = (5, 6, 7)

        # energy bounds matching those of the
        # generated weight windows
        e_bnds = [0.0, 0.5, 2E7]

        ww_n = openmc.WeightWindows(ww_mesh,
                                    ww_n_lower_bnds,
                                    None,
                                    10.0,
                                    e_bnds,
                                    survival_ratio=1.01)

        ww_p = openmc.WeightWindows(ww_mesh,
                                    ww_p_lower_bnds,
                                    None,
                                    10.0,
                                    e_bnds,
                                    survival_ratio=1.01)

        model.settings.weight_windows = [ww_n, ww_p]

        # check that string form of the class can be created
        for ww in model.settings.weight_windows:
            str(ww)

        # run again with variance reduction on
        model.settings.weight_windows_on = True
        ww_sp = model.run()
        os.rename(ww_sp, 'statepoint.ww.h5')

        # load both statepoints and examine results
        asp = openmc.StatePoint('statepoint.analog.h5')
        wsp = openmc.StatePoint('statepoint.ww.h5')

        analog_tally = asp.tallies[1]
        ww_tally = wsp.tallies[1]

        def compare_results(particle, analog_tally, ww_tally):
            # get values from each of the tallies
            an_mean = analog_tally.get_values(filters=[openmc.ParticleFilter],
                                              filter_bins=[(particle,)])
            ww_mean = ww_tally.get_values(filters=[openmc.ParticleFilter],
                                          filter_bins=[(particle,)])

            # expect that more bins were scored with weight windows than
            # the analog run
            assert np.count_nonzero(an_mean) < np.count_nonzero(ww_mean)

            an_rel_err = analog_tally.get_values(filters=[openmc.ParticleFilter],
                                                 filter_bins=[(particle,)],
                                                 value='rel_err')
            ww_rel_err = ww_tally.get_values(filters=[openmc.ParticleFilter],
                                             filter_bins=[(particle,)],
                                             value='rel_err')

            an_rel_err[an_mean == 0.0] = 1.0
            ww_rel_err[ww_mean == 0.0] = 1.0

            an_avg_rel_err = np.mean(an_rel_err)
            ww_avg_rel_err = np.mean(ww_rel_err)

            # expect that the average relative error in the tally
            # decreases
            assert an_avg_rel_err > ww_avg_rel_err

            # ensure that the value of the mesh bin containing the
            # source is statistically similar in both runs
            an_std_dev = analog_tally.get_values(filters=[openmc.ParticleFilter],
                                                 filter_bins=[(particle,)],
                                                 value='std_dev')
            ww_std_dev = ww_tally.get_values(filters=[openmc.ParticleFilter],
                                             filter_bins=[(particle,)],
                                             value='std_dev')

            # index of the mesh bin containing the source for the higher
            # energy group
            source_bin_idx = (an_mean.shape[0]//2, 0, 0)

            an_source_bin = ufloat(an_mean[source_bin_idx],
                                   an_std_dev[source_bin_idx])
            ww_source_bin = ufloat(ww_mean[source_bin_idx],
                                   ww_std_dev[source_bin_idx])

            diff = an_source_bin - ww_source_bin

            # check that values are within two combined standard deviations
            assert abs(diff.nominal_value) / diff.std_dev < 2.0

        compare_results('neutron', analog_tally, ww_tally)
        compare_results('photon', analog_tally, ww_tally)
