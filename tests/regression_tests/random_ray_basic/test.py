import os

import numpy as np
import openmc

from tests.testing_harness import TolerantPyAPITestHarness


class MGXSTestHarness(TolerantPyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def random_ray_model() -> openmc.Model:
    ###############################################################################
    # Create multigroup data

    # Instantiate the energy group data
    group_edges = [1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges)

    # Instantiate the 7-group (C5G7) cross section data
    uo2_xsdata = openmc.XSdata('UO2', groups)
    uo2_xsdata.order = 0
    uo2_xsdata.set_total(
        [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
         0.5644058])
    uo2_xsdata.set_absorption([8.0248e-03, 3.7174e-03, 2.6769e-02, 9.6236e-02,
                               3.0020e-02, 1.1126e-01, 2.8278e-01])
    scatter_matrix = np.array(
        [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    uo2_xsdata.set_scatter_matrix(scatter_matrix)
    uo2_xsdata.set_fission([7.21206e-03, 8.19301e-04, 6.45320e-03,
                            1.85648e-02, 1.78084e-02, 8.30348e-02,
                            2.16004e-01])
    uo2_xsdata.set_nu_fission([2.005998e-02, 2.027303e-03, 1.570599e-02,
                               4.518301e-02, 4.334208e-02, 2.020901e-01,
                               5.257105e-01])
    uo2_xsdata.set_chi([5.8791e-01, 4.1176e-01, 3.3906e-04, 1.1761e-07, 0.0000e+00,
                        0.0000e+00, 0.0000e+00])

    h2o_xsdata = openmc.XSdata('LWTR', groups)
    h2o_xsdata.order = 0
    h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
                          0.718, 1.2544497, 2.650379])
    h2o_xsdata.set_absorption([6.0105e-04, 1.5793e-05, 3.3716e-04,
                               1.9406e-03, 5.7416e-03, 1.5001e-02,
                               3.7239e-02])
    scatter_matrix = np.array(
        [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
          [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
          [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
          [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
          [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
          [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    h2o_xsdata.set_scatter_matrix(scatter_matrix)

    mg_cross_sections = openmc.MGXSLibrary(groups)
    mg_cross_sections.add_xsdatas([uo2_xsdata, h2o_xsdata])
    mg_cross_sections.export_to_hdf5()

    ###############################################################################
    # Create materials for the problem

    # Instantiate some Materials and register the appropriate macroscopic data
    uo2 = openmc.Material(name='UO2 fuel')
    uo2.set_density('macro', 1.0)
    uo2.add_macroscopic('UO2')

    water = openmc.Material(name='Water')
    water.set_density('macro', 1.0)
    water.add_macroscopic('LWTR')

    # Instantiate a Materials collection and export to XML
    materials = openmc.Materials([uo2, water])
    materials.cross_sections = "mgxs.h5"

    ###############################################################################
    # Define problem geometry

    ########################################
    # Define an unbounded pincell universe

    pitch = 1.26

    # Create a surface for the fuel outer radius
    fuel_or = openmc.ZCylinder(r=0.54, name='Fuel OR')
    inner_ring_a = openmc.ZCylinder(r=0.33, name='inner ring a')
    inner_ring_b = openmc.ZCylinder(r=0.45, name='inner ring b')
    outer_ring_a = openmc.ZCylinder(r=0.60, name='outer ring a')
    outer_ring_b = openmc.ZCylinder(r=0.69, name='outer ring b')

    # Instantiate Cells
    fuel_a = openmc.Cell(fill=uo2, region=-inner_ring_a, name='fuel inner a')
    fuel_b = openmc.Cell(fill=uo2, region=+inner_ring_a & -inner_ring_b, name='fuel inner b')
    fuel_c = openmc.Cell(fill=uo2, region=+inner_ring_b & -fuel_or, name='fuel inner c')
    moderator_a = openmc.Cell(fill=water, region=+fuel_or & -outer_ring_a, name='moderator inner a')
    moderator_b = openmc.Cell(fill=water, region=+outer_ring_a & -outer_ring_b, name='moderator outer b')
    moderator_c = openmc.Cell(fill=water, region=+outer_ring_b, name='moderator outer c')

    # Create pincell universe
    pincell_base = openmc.Universe()

    # Register Cells with Universe
    pincell_base.add_cells([fuel_a, fuel_b, fuel_c, moderator_a, moderator_b, moderator_c])

    # Create planes for azimuthal sectors
    azimuthal_planes = []
    for i in range(8):
        angle = 2 * i * openmc.pi / 8
        normal_vector = (-openmc.sin(angle), openmc.cos(angle), 0)
        azimuthal_planes.append(openmc.Plane(a=normal_vector[0], b=normal_vector[1], c=normal_vector[2], d=0))

    # Create a cell for each azimuthal sector
    azimuthal_cells = []
    for i in range(8):
        azimuthal_cell = openmc.Cell(name=f'azimuthal_cell_{i}')
        azimuthal_cell.fill = pincell_base
        azimuthal_cell.region = +azimuthal_planes[i] & -azimuthal_planes[(i+1) % 8]
        azimuthal_cells.append(azimuthal_cell)

    # Create a geometry with the azimuthal universes
    pincell = openmc.Universe(cells=azimuthal_cells)

    ########################################
    # Define a moderator lattice universe

    moderator_infinite = openmc.Cell(fill=water, name='moderator infinite')
    mu = openmc.Universe()
    mu.add_cells([moderator_infinite])

    lattice = openmc.RectLattice()
    lattice.lower_left = [-pitch/2.0, -pitch/2.0]
    lattice.pitch = [pitch/10.0, pitch/10.0]
    lattice.universes = np.full((10, 10), mu)

    mod_lattice_cell = openmc.Cell(fill=lattice)

    mod_lattice_uni = openmc.Universe()

    mod_lattice_uni.add_cells([mod_lattice_cell])

    ########################################
    # Define 2x2 outer lattice
    lattice2x2 = openmc.RectLattice()
    lattice2x2.lower_left = (-pitch, -pitch)
    lattice2x2.pitch = (pitch, pitch)
    lattice2x2.universes = [
        [pincell, pincell],
        [pincell, mod_lattice_uni]
    ]

    ########################################
    # Define cell containing lattice and other stuff
    box = openmc.model.RectangularPrism(pitch*2, pitch*2, boundary_type='reflective')

    assembly = openmc.Cell(fill=lattice2x2, region=-box, name='assembly')

    # Create a geometry with the top-level cell
    geometry = openmc.Geometry([assembly])

    ###############################################################################
    # Define problem settings

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.batches = 10
    settings.inactive = 5
    settings.particles = 100

    # Create an initial uniform spatial source distribution over fissionable zones
    lower_left = (-pitch, -pitch, -1)
    upper_right = (pitch, pitch, 1)
    uniform_dist = openmc.stats.Box(lower_left, upper_right)
    rr_source = openmc.IndependentSource(space=uniform_dist)

    settings.random_ray['distance_active'] = 100.0
    settings.random_ray['distance_inactive'] = 20.0
    settings.random_ray['ray_source'] = rr_source
    settings.random_ray['volume_normalized_flux_tallies'] = True

    ###############################################################################
    # Define tallies

    # Create a mesh that will be used for tallying
    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 2)
    mesh.lower_left = (-pitch, -pitch)
    mesh.upper_right = (pitch, pitch)

    # Create a mesh filter that can be used in a tally
    mesh_filter = openmc.MeshFilter(mesh)

    # Create an energy group filter as well
    energy_filter = openmc.EnergyFilter(group_edges)

    # Now use the mesh filter in a tally and indicate what scores are desired
    tally = openmc.Tally(name="Mesh tally")
    tally.filters = [mesh_filter, energy_filter]
    tally.scores = ['flux', 'fission', 'nu-fission']
    tally.estimator = 'analog'

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([tally])

    ###############################################################################
    #                   Exporting to OpenMC model
    ###############################################################################

    model = openmc.Model()
    model.geometry = geometry
    model.materials = materials
    model.settings = settings
    model.tallies = tallies
    return model


def test_random_ray_basic():
    harness = MGXSTestHarness('statepoint.10.h5', random_ray_model())
    harness.main()
