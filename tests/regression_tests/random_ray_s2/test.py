import os

import openmc
import openmc.model
import numpy as np

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_random_ray_s2():
    NUM_SOURCE_REGIONS = 10
    L = 40.0

    # Make the simple MGXS library and material.
    groups = openmc.mgxs.EnergyGroups(group_edges=[1e-5, 20.0e6])

    domain_mat_data = openmc.XSdata('domain', groups)
    domain_mat_data.order = 0
    domain_mat_data.set_total([0.1])
    domain_mat_data.set_absorption([0.1])
    domain_mat_data.set_scatter_matrix(np.rollaxis(np.array([[[0.0]]]), 0, 3))
    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas([domain_mat_data])
    mg_cross_sections_file.export_to_hdf5()

    domain_data = openmc.Macroscopic('domain')
    domain_mat = openmc.Material(name='domain')
    domain_mat.set_density('macro', 1.0)
    domain_mat.add_macroscopic(domain_data)

    container = openmc.model.RectangularPrism(width=10.0, height=10.0, axis='x',
                                              origin=(0.0, 0.0), boundary_type='reflective')
    left = openmc.XPlane(x0 = 0.0, boundary_type='vacuum')
    right = openmc.XPlane(x0 = L, boundary_type='vacuum')
    cell = [openmc.Cell(region = +left & -right & -container, fill = domain_mat)]

    model = openmc.model.Model()
    model.geometry = openmc.Geometry(root=openmc.Universe(cells=cell))
    model.materials = openmc.Materials([domain_mat])
    model.materials.cross_sections = './mgxs.h5'

    mesh = openmc.RegularMesh()
    mesh.dimension = (NUM_SOURCE_REGIONS, 1, 1)
    mesh.lower_left = (0.0, -5.0, -5.0)
    mesh.upper_right = (L, 5.0, 5.0)

    tally = openmc.Tally(name="LR")
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ['flux']
    tally.estimator = 'tracklength'
    model.tallies.append(tally)

    uniform_dist = openmc.stats.Box((0.0, -5.0, -5.0), (L, 5.0, 5.0))
    model.settings.source = [
      openmc.IndependentSource(space=uniform_dist,
                              energy=openmc.stats.Discrete(x = 1e3, p = 1.0),
                              constraints={'domains' : [domain_mat]})
    ]
    model.settings.energy_mode = "multi-group"
    model.settings.batches = 30
    model.settings.inactive = 10
    model.settings.particles = 100
    model.settings.run_mode = 'fixed source'
    model.settings.random_ray['distance_inactive'] = 100.0
    model.settings.random_ray['distance_active'] = 400.0
    model.settings.random_ray['ray_source'] = openmc.IndependentSource(space=uniform_dist)
    model.settings.random_ray['source_shape'] = 'flat'
    model.settings.random_ray['sample_method'] = 's2'
    model.settings.random_ray['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]

    model.export_to_model_xml()

    harness = MGXSTestHarness('statepoint.30.h5', model)
    harness.main()
