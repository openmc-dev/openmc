import os

import numpy as np
import openmc
from openmc.examples import random_ray_three_region_cube

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_low_density():
    model = random_ray_three_region_cube()

    # Rebuild the MGXS library to have a material with very
    # low macroscopic cross sections
    ebins = [1e-5, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

    void_sigma_a = 4.0e-6
    void_sigma_s = 3.0e-4
    void_mat_data = openmc.XSdata('void', groups)
    void_mat_data.order = 0
    void_mat_data.set_total([void_sigma_a + void_sigma_s])
    void_mat_data.set_absorption([void_sigma_a])
    void_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[void_sigma_s]]]), 0, 3))

    absorber_sigma_a = 0.75
    absorber_sigma_s = 0.25
    absorber_mat_data = openmc.XSdata('absorber', groups)
    absorber_mat_data.order = 0
    absorber_mat_data.set_total([absorber_sigma_a + absorber_sigma_s])
    absorber_mat_data.set_absorption([absorber_sigma_a])
    absorber_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[absorber_sigma_s]]]), 0, 3))

    multiplier = 0.0000001
    source_sigma_a = void_sigma_a * multiplier
    source_sigma_s = void_sigma_s * multiplier
    source_mat_data = openmc.XSdata('source', groups)
    source_mat_data.order = 0
    source_mat_data.set_total([source_sigma_a + source_sigma_s])
    source_mat_data.set_absorption([source_sigma_a])
    source_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[source_sigma_s]]]), 0, 3))

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas(
        [source_mat_data, void_mat_data, absorber_mat_data])
    mg_cross_sections_file.export_to_hdf5()

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
