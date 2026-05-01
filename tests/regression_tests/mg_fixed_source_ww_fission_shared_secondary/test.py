import os

import numpy as np
import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    # Instantiate the energy group data and file object
    groups = openmc.mgxs.EnergyGroups([0.0, 0.625, 20.0e6])

    mg_cross_sections_file = openmc.MGXSLibrary(groups)

    # Make the base, isotropic data
    nu = [2.50, 2.50]
    fiss = np.array([0.002817, 0.097])
    capture = [0.008708, 0.02518]
    absorption = np.add(capture, fiss)
    scatter = np.array(
        [[[0.31980, 0.06694], [0.004555, -0.0003972]],
         [[0.00000, 0.00000], [0.424100, 0.05439000]]])
    total = [0.33588, 0.54628]
    chi = [1., 0.]

    mat_1 = openmc.XSdata('mat_1', groups)
    mat_1.order = 1
    mat_1.set_nu_fission(np.multiply(nu, fiss))
    mat_1.set_absorption(absorption)
    mat_1.set_scatter_matrix(scatter)
    mat_1.set_total(total)
    mat_1.set_chi(chi)
    mg_cross_sections_file.add_xsdata(mat_1)

    # Write the file
    mg_cross_sections_file.export_to_hdf5('2g.h5')


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_fixed_source_ww_fission_shared_secondary():
    create_library()
    model = slab_mg()

    # Override settings for fixed-source mode with shared secondary bank
    model.settings.run_mode = 'fixed source'
    model.settings.inactive = 0
    model.settings.batches = 2
    model.settings.particles = 100
    model.settings.create_fission_neutrons = True
    model.settings.shared_secondary_bank = True
    model.settings.max_history_splits = 100

    # Add weight windows on a simple 1D mesh
    ww_mesh = openmc.RegularMesh()
    ww_mesh.lower_left = (0.0, -1000.0, -1000.0)
    ww_mesh.upper_right = (929.45, 1000.0, 1000.0)
    ww_mesh.dimension = (5, 1, 1)

    # Uniform lower bounds for 2 energy groups, 5 spatial bins
    lower_bounds = np.full((2, 5, 1, 1), 0.5)
    ww = openmc.WeightWindows(
        ww_mesh,
        lower_bounds.flatten(),
        None,
        5.0,
        [0.0, 0.625, 20.0e6],
        'neutron'
    )
    model.settings.weight_windows = [ww]

    # Add a flux tally
    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, -1000.0, -1000.0)
    mesh.upper_right = (929.45, 1000.0, 1000.0)
    mesh.dimension = (5, 1, 1)

    tally = openmc.Tally()
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ['flux']
    model.tallies = [tally]

    harness = MGXSTestHarness('statepoint.2.h5', model)
    harness.main()
