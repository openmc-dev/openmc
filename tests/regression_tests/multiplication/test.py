import openmc
from openmc.stats import delta_function
import numpy as np
import pytest
from openmc.examples import slab_mg
import os
    
from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)

@pytest.fixture()
def slab_model():
    model = slab_mg(mgxslib_name='mgxs.h5')
    right_boundary = model.geometry.get_all_surfaces()[2]
    right_boundary.coefficients['x0'] = 10.0
    
    cell = model.geometry.get_all_cells()[1]
    
    mat = model.geometry.get_all_materials()[1]
    mat.set_density('macro', 0.01) 
    ###########################################################################
    # Create multigroup data

    # Instantiate the energy group data
    ebins = np.geomspace(1e-5, 20.0e6, 5)
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

    nusigma_f = np.array([9.6,5.4,5.2,2.5])
    sigma_s = np.array([[0.5,0.5,0.5,0.5],
                        [0.0,1.0,0.5,0.5],
                        [0.0,0.0,1.5,0.5],
                        [0.0,0.0,0.0,2.0]])
    sigma_t = np.array([5.0,5.0,5.0,5.0])
    sigma_a = sigma_t - sigma_s.sum(axis=1)
    chi = np.array([0.0,0.2,0.8,0.0])
    mat_data = openmc.XSdata('mat_1', groups)
    mat_data.order = 0
    mat_data.set_total(sigma_t)
    mat_data.set_nu_fission(nusigma_f)
    mat_data.set_chi(chi)
    mat_data.set_absorption(sigma_a)
    mat_data.set_scatter_matrix(sigma_s[...,np.newaxis])

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdata(mat_data)
    mg_cross_sections_file.export_to_hdf5()

    # Settings
    model.settings.particles = 100000
    model.settings.batches = 20
    model.settings.inactive = 10

    space = openmc.stats.Box([0,-1000,-1000],[10,1000,1000])
    model.settings.source = openmc.IndependentSource(
        space=space, energy = delta_function(10e6))
        
    return model


def test_multiplication(slab_model):
    harness = MGXSTestHarness("statepoint.20.h5", model=slab_model)
    harness.main()
