import os

import numpy as np

import openmc

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


def create_model():
    create_library()

    # # Make Materials
    materials_file = openmc.Materials()

    mat_names = ['base leg', 'base tab', 'base hist', 'base matrix', 'base ang']
    macros = []
    mats = []
    for i in range(len(mat_names)):
        macros.append(openmc.Macroscopic('mat_' + str(i + 1)))
        mats.append(openmc.Material(name=mat_names[i]))
        mats[-1].set_density('macro', 1.0)
        mats[-1].add_macroscopic(macros[-1])

    # Add in the microscopic data
    mats.append(openmc.Material(name='micro'))
    mats[-1].set_density("sum")
    mats[-1].add_nuclide("mat_1", 0.5)
    mats[-1].add_nuclide("mat_6", 0.5)

    materials_file += mats

    materials_file.cross_sections = '2g.h5'

    # # Make Geometry
    rad_outer = 929.45
    # Set a cell boundary to exist for every material above (exclude the 0)
    rads = np.linspace(0., rad_outer, len(mats) + 1, endpoint=True)[1:]

    # Instantiate Universe
    root = openmc.Universe(universe_id=0, name='root universe')
    cells = []

    surfs = []
    surfs.append(openmc.XPlane(x0=0., boundary_type='reflective'))
    for r, rad in enumerate(rads):
        if r == len(rads) - 1:
            surfs.append(openmc.XPlane(x0=rad, boundary_type='vacuum'))
        else:
            surfs.append(openmc.XPlane(x0=rad))

    # Instantiate Cells
    cells = []
    for c in range(len(surfs) - 1):
        cells.append(openmc.Cell())
        cells[-1].region = (+surfs[c] & -surfs[c + 1])
        cells[-1].fill = mats[c]

    # Register Cells with Universe
    root.add_cells(cells)

    # Instantiate a Geometry, register the root Universe, and export to XML
    geometry_file = openmc.Geometry(root)

    # # Make Settings
    # Instantiate a Settings object, set all runtime parameters
    settings_file = openmc.Settings()
    settings_file.energy_mode = "multi-group"
    settings_file.tabular_legendre = {'enable': False}
    settings_file.batches = 10
    settings_file.inactive = 5
    settings_file.particles = 1000

    # Build source distribution
    INF = 1000.
    bounds = [0., -INF, -INF, rads[0], INF, INF]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings_file.source = openmc.source.Source(space=uniform_dist)

    settings_file.output = {'summary': False}

    model = openmc.model.Model()
    model.geometry = geometry_file
    model.materials = materials_file
    model.settings = settings_file

    return model


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = '2g .h5'
        if os.path.exists(f):
            os.remove(f)


def test_mg_benchmark():
    model = create_model()
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
