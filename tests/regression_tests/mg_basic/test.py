import os

import numpy as np

import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    # Instantiate the energy group data and file object
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])

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

    # Make a version of mat-1 which has a tabular representation of the
    # scattering vice Legendre with 33 points
    mat_2 = mat_1.convert_scatter_format('tabular', 33)
    mat_2.name = 'mat_2'
    mg_cross_sections_file.add_xsdata(mat_2)

    # Make a version of mat-1 which has a histogram representation of the
    # scattering vice Legendre with 33 bins
    mat_3 = mat_1.convert_scatter_format('histogram', 33)
    mat_3.name = 'mat_3'
    mg_cross_sections_file.add_xsdata(mat_3)

    # Make a version which uses a fission matrix vice chi & nu-fission
    mat_4 = openmc.XSdata('mat_4', groups)
    mat_4.order = 1
    mat_4.set_nu_fission(np.outer(np.multiply(nu, fiss), chi))
    mat_4.set_absorption(absorption)
    mat_4.set_scatter_matrix(scatter)
    mat_4.set_total(total)
    mg_cross_sections_file.add_xsdata(mat_4)

    # Make an angle-dependent version of mat_1 with 2 polar and 2 azim. angles
    mat_5 = mat_1.convert_representation('angle', 2, 2)
    mat_5.name = 'mat_5'
    mg_cross_sections_file.add_xsdata(mat_5)

    # Make a copy of mat_1 for testing microscopic cross sections
    mat_6 = openmc.XSdata('mat_6', groups)
    mat_6.order = 1
    mat_6.set_nu_fission(np.multiply(nu, fiss))
    mat_6.set_absorption(absorption)
    mat_6.set_scatter_matrix(scatter)
    mat_6.set_total(total)
    mat_6.set_chi(chi)
    mg_cross_sections_file.add_xsdata(mat_6)

    # Write the file
    mg_cross_sections_file.export_to_hdf5('2g.h5')


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_basic():
    create_library()
    mat_names = ['base leg', 'base tab', 'base hist', 'base matrix',
                 'base ang', 'micro']
    model = slab_mg(num_regions=6, mat_names=mat_names)
    # Modify the last material to be a microscopic combination of nuclides
    model.materials[-1] = openmc.Material(name='micro', material_id=6)
    model.materials[-1].set_density("sum")
    model.materials[-1].add_nuclide("mat_1", 0.5)
    model.materials[-1].add_nuclide("mat_6", 0.5)

    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
