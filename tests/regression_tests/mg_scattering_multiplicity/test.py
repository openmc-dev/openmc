import os

import numpy as np

import openmc
from openmc.examples import slab_mg

from tests.testing_harness import PyAPITestHarness


def create_library():
    groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625, 20.0e6])

    # Two-group data in which the scattering multiplicity depends on the
    # outgoing group: every fast-to-thermal transfer produces two neutrons
    scatter = np.array([[0.2, 0.2], [0.0, 0.2]])
    multiplicity = np.array([[1.0, 2.0], [1.0, 1.0]])

    mat_1 = openmc.XSdata('mat_1', groups)
    mat_1.order = 0
    mat_1.set_total([1.0, 1.0])
    mat_1.set_absorption([0.6, 0.8])
    mat_1.set_fission([0.16, 0.32])
    mat_1.set_nu_fission([0.4, 0.8])
    mat_1.set_chi([1.0, 0.0])
    # With multiplicity data, the scattering matrix holds neutron production
    mat_1.set_scatter_matrix((scatter * multiplicity)[:, :, np.newaxis])
    mat_1.set_multiplicity_matrix(multiplicity)

    # Tabular and histogram representations of the same data, so that each of
    # the three scattering samplers is exercised
    mat_2 = mat_1.convert_scatter_format('tabular', 33)
    mat_2.name = 'mat_2'
    mat_3 = mat_1.convert_scatter_format('histogram', 2)
    mat_3.name = 'mat_3'

    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    for mat in (mat_1, mat_2, mat_3):
        mg_cross_sections_file.add_xsdata(mat)
    mg_cross_sections_file.export_to_hdf5('2g.h5')


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_scattering_multiplicity():
    create_library()
    mat_names = ['legendre', 'tabular', 'histogram']
    model = slab_mg(num_regions=3, mat_names=mat_names)
    model.settings.particles = 100

    # Start particles in all three regions (the default source of slab_mg only
    # covers the first one)
    x_max = model.geometry.bounding_box.upper_right[0]
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
        [0.0, -1000.0, -1000.0], [x_max, 1000.0, 1000.0]))

    # The scatter score uses the total scattering-event cross section, which
    # also depends on the multiplicity data
    tally = openmc.Tally()
    tally.filters = [openmc.MaterialFilter(model.materials),
                     openmc.EnergyFilter([0.0, 0.625, 20.0e6])]
    tally.scores = ['flux', 'scatter', 'nu-scatter']
    model.tallies = [tally]

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
