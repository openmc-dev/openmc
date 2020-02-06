import os

import numpy as np

import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    # Instantiate the energy group data and file object
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])

    mg_cross_sections_file = openmc.MGXSLibrary(groups, 6)

    # Make the base, isotropic data
    nu = np.array([2.50, 2.50])
    fiss = np.array([0.002817, 0.097])
    capture = np.array([0.008708, 0.02518])
    absorption = capture + fiss
    scatter = np.array(
        [[[0.31980, 0.06694], [0.004555, -0.0003972]],
         [[0.00000, 0.00000], [0.424100, 0.05439000]]])
    total = np.array([0.33588, 0.54628])
    chi = np.array([1., 0.])
    decay_rate = np.array([0.013336, 0.032739, 0.12078, 0.30278, 0.84949,
                           2.853])
    delayed_yield = np.array([0.00055487, 0.00286407, 0.00273429, 0.0061305,
                              0.00251342, 0.00105286])
    inv_vel = 1.0 / np.array([1.4e9, 4.4e5])


    mat_1 = openmc.XSdata('mat_1', groups, num_delayed_groups=6)
    mat_1.order = 1
    mat_1.set_fission(fiss)
    mat_1.set_kappa_fission(fiss * 200e6)
    mat_1.set_nu_fission(nu * fiss)
    mat_1.set_beta(delayed_yield / 2.5)
    mat_1.set_decay_rate(decay_rate)
    mat_1.set_absorption(absorption)
    mat_1.set_scatter_matrix(scatter)
    mat_1.set_total(total)
    mat_1.set_chi(chi)
    mat_1.set_inverse_velocity(inv_vel)
    mg_cross_sections_file.add_xsdata(mat_1)

    # Write the file
    mg_cross_sections_file.export_to_hdf5('2g.h5')


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_tallies():
    create_library()
    model = slab_mg()

    # Instantiate a tally mesh
    mesh = openmc.RegularMesh(mesh_id=1)
    mesh.dimension = [10, 1, 1]
    mesh.lower_left = [0.0, 0.0, 0.0]
    mesh.upper_right = [929.45, 1000, 1000]

    # Instantiate some tally filters
    energy_filter = openmc.EnergyFilter([0.0, 20.0e6])
    energyout_filter = openmc.EnergyoutFilter([0.0, 20.0e6])
    energies = [0.0, 0.625, 20.0e6]
    matching_energy_filter = openmc.EnergyFilter(energies)
    matching_eout_filter = openmc.EnergyoutFilter(energies)
    mesh_filter = openmc.MeshFilter(mesh)

    mat_filter = openmc.MaterialFilter(model.materials)

    nuclides = model.xs_data

    scores_with_nuclides = [
      'total', 'absorption', 'fission', 'nu-fission', 'inverse-velocity',
      'prompt-nu-fission', 'delayed-nu-fission', 'kappa-fission', 'events',
      'decay-rate']
    scores_without_nuclides = scores_with_nuclides + ['flux']

    for do_nuclides, scores in ((False, scores_without_nuclides),
                                (True, scores_with_nuclides)):
        t = openmc.Tally()
        t.filters = [mesh_filter]
        t.estimator = 'analog'
        t.scores = scores
        if do_nuclides:
            t.nuclides = nuclides
        model.tallies.append(t)

        t = openmc.Tally()
        t.filters = [mesh_filter]
        t.estimator = 'tracklength'
        t.scores = scores
        if do_nuclides:
            t.nuclides = nuclides
        model.tallies.append(t)

        # Impose energy bins that dont match the MG structure and those
        # that do
        for match_energy_bins in [False, True]:
            if match_energy_bins:
                e_filter = matching_energy_filter
                eout_filter = matching_eout_filter
            else:
                e_filter = energy_filter
                eout_filter = energyout_filter

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'analog'
            t.scores = scores + ['scatter', 'nu-scatter']
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'collision'
            t.scores = scores
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter]
            t.estimator = 'tracklength'
            t.scores = scores
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

            t = openmc.Tally()
            t.filters = [mat_filter, e_filter, eout_filter]
            t.scores = ['scatter', 'nu-scatter', 'nu-fission']
            if do_nuclides:
                t.nuclides = nuclides
            model.tallies.append(t)

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
