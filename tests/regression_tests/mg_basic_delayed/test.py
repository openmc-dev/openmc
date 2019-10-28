import os

import numpy as np

import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    # Instantiate the energy group data and file object
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])
    n_dg = 2

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.num_delayed_groups = n_dg

    beta = np.array([0.003, 0.003])
    one_m_beta = 1. - np.sum(beta)
    nu = [2.50, 2.50]
    fiss = np.array([0.002817, 0.097])
    capture = [0.008708, 0.02518]
    absorption = np.add(capture, fiss)
    scatter = np.array(
        [[[0.31980, 0.06694], [0.004555, -0.0003972]],
         [[0.00000, 0.00000], [0.424100, 0.05439000]]])
    total = [0.33588, 0.54628]
    chi = [1., 0.]

    # Make the base data that uses chi & nu-fission vectors with a beta
    mat_1 = openmc.XSdata('mat_1', groups)
    mat_1.order = 1
    mat_1.num_delayed_groups = 2
    mat_1.set_beta(beta)
    mat_1.set_nu_fission(np.multiply(nu, fiss))
    mat_1.set_absorption(absorption)
    mat_1.set_scatter_matrix(scatter)
    mat_1.set_total(total)
    mat_1.set_chi(chi)
    mg_cross_sections_file.add_xsdata(mat_1)

    # Make a version that uses prompt and delayed version of nufiss and chi
    mat_2 = openmc.XSdata('mat_2', groups)
    mat_2.order = 1
    mat_2.num_delayed_groups = 2
    mat_2.set_prompt_nu_fission(one_m_beta * np.multiply(nu, fiss))
    delay_nu_fiss = np.zeros((n_dg, groups.num_groups))
    for dg in range(n_dg):
        for g in range(groups.num_groups):
            delay_nu_fiss[dg, g] = beta[dg] * nu[g] * fiss[g]
    mat_2.set_delayed_nu_fission(delay_nu_fiss)
    mat_2.set_absorption(absorption)
    mat_2.set_scatter_matrix(scatter)
    mat_2.set_total(total)
    mat_2.set_chi_prompt(chi)
    mat_2.set_chi_delayed(np.stack([chi] * n_dg))
    mg_cross_sections_file.add_xsdata(mat_2)

    # Make a version that uses a nu-fission matrix with a beta
    mat_3 = openmc.XSdata('mat_3', groups)
    mat_3.order = 1
    mat_3.num_delayed_groups = 2
    mat_3.set_beta(beta)
    mat_3.set_nu_fission(np.outer(np.multiply(nu, fiss), chi))
    mat_3.set_absorption(absorption)
    mat_3.set_scatter_matrix(scatter)
    mat_3.set_total(total)
    mg_cross_sections_file.add_xsdata(mat_3)

    # Make a version that uses prompt and delayed version of the nufiss matrix
    mat_4 = openmc.XSdata('mat_4', groups)
    mat_4.order = 1
    mat_4.num_delayed_groups = 2
    mat_4.set_prompt_nu_fission(one_m_beta *
                                np.outer(np.multiply(nu, fiss), chi))
    delay_nu_fiss = np.zeros((n_dg, groups.num_groups, groups.num_groups))
    for dg in range(n_dg):
        for g in range(groups.num_groups):
            for go in range(groups.num_groups):
                delay_nu_fiss[dg, g, go] = beta[dg] * nu[g] * fiss[g] * chi[go]
    mat_4.set_delayed_nu_fission(delay_nu_fiss)
    mat_4.set_absorption(absorption)
    mat_4.set_scatter_matrix(scatter)
    mat_4.set_total(total)
    mg_cross_sections_file.add_xsdata(mat_4)

    # Make the base data that uses chi & nu-fiss vectors with a group-wise beta
    mat_5 = openmc.XSdata('mat_5', groups)
    mat_5.order = 1
    mat_5.num_delayed_groups = 2
    mat_5.set_beta(np.stack([beta] * groups.num_groups))
    mat_5.set_nu_fission(np.multiply(nu, fiss))
    mat_5.set_absorption(absorption)
    mat_5.set_scatter_matrix(scatter)
    mat_5.set_total(total)
    mat_5.set_chi(chi)
    mg_cross_sections_file.add_xsdata(mat_5)

    # Make a version that uses a nu-fission matrix with a group-wise beta
    mat_6 = openmc.XSdata('mat_6', groups)
    mat_6.order = 1
    mat_6.num_delayed_groups = 2
    mat_6.set_beta(np.stack([beta] * groups.num_groups))
    mat_6.set_nu_fission(np.outer(np.multiply(nu, fiss), chi))
    mat_6.set_absorption(absorption)
    mat_6.set_scatter_matrix(scatter)
    mat_6.set_total(total)
    mg_cross_sections_file.add_xsdata(mat_6)

    # Write the file
    mg_cross_sections_file.export_to_hdf5('2g.h5')


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_basic_delayed():
    create_library()
    model = slab_mg(num_regions=6, mat_names=['vec beta', 'vec no beta',
                                              'matrix beta', 'matrix no beta',
                                              'vec group beta',
                                              'matrix group beta'])

    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
