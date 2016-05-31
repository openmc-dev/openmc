import openmc
import openmc.mgxs

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 100
inactive = 10
particles = 1000

###############################################################################
#                 Exporting to OpenMC mgxs.xml file
###############################################################################

# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(group_edges=[1E-11, 0.0635E-6, 10.0E-6,
                                               1.0E-4, 1.0E-3, 0.5, 1.0, 20.0])

# Instantiate the 7-group (C5G7) cross section data
uo2_xsdata = openmc.XSdata('UO2.300K', groups)
uo2_xsdata.order = 0
uo2_xsdata.total = [0.1779492, 0.3298048, 0.4803882, 0.5543674,
                    0.3118013, 0.3951678, 0.5644058]
uo2_xsdata.absorption = [8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02,
                         3.0020E-02, 1.1126E-01, 2.8278E-01]
uo2_xsdata.scatter = [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
                       [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
                       [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
                       [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
                       [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
                       [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
                       [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]]
uo2_xsdata.fission = [7.21206E-03, 8.19301E-04, 6.45320E-03,
                      1.85648E-02, 1.78084E-02, 8.30348E-02,
                      2.16004E-01]
uo2_xsdata.nu_fission = [2.005998E-02, 2.027303E-03, 1.570599E-02,
                         4.518301E-02, 4.334208E-02, 2.020901E-01,
                         5.257105E-01]
uo2_xsdata.chi = [5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07,
                  0.0000E+00, 0.0000E+00, 0.0000E+00]

h2o_xsdata = openmc.XSdata('LWTR.300K', groups)
h2o_xsdata.order = 0
h2o_xsdata.total = [0.15920605, 0.412969593, 0.59030986, 0.58435,
                    0.718, 1.2544497, 2.650379]
h2o_xsdata.absorption = [6.0105E-04, 1.5793E-05, 3.3716E-04,
                         1.9406E-03, 5.7416E-03, 1.5001E-02,
                         3.7239E-02]
h2o_xsdata.scatter = [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
                       [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
                       [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
                       [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
                       [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
                       [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
                       [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]]

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
mg_cross_sections_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate some Macroscopic Data
uo2_data = openmc.Macroscopic('UO2', '300K')
h2o_data = openmc.Macroscopic('LWTR', '300K')

# Instantiate some Materials and register the appropriate Macroscopic objects
uo2 = openmc.Material(material_id=1, name='UO2 fuel')
uo2.set_density('macro', 1.0)
uo2.add_macroscopic(uo2_data)

water = openmc.Material(material_id=2, name='Water')
water.set_density('macro', 1.0)
water.add_macroscopic(h2o_data)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, water])
materials_file.default_xs = '300K'
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate ZCylinder surfaces
fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=0.54, name='Fuel OR')
left = openmc.XPlane(surface_id=4, x0=-0.63, name='left')
right = openmc.XPlane(surface_id=5, x0=0.63, name='right')
bottom = openmc.YPlane(surface_id=6, y0=-0.63, name='bottom')
top = openmc.YPlane(surface_id=7, y0=0.63, name='top')

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
top.boundary_type = 'reflective'
bottom.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(cell_id=1, name='cell 1')
moderator = openmc.Cell(cell_id=2, name='cell 2')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
moderator.region = +fuel_or & +left & -right & +bottom & -top

# Register Materials with Cells
fuel.fill = uo2
moderator.fill = water

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
root.add_cells([fuel, moderator])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.energy_mode = "multi-group"
settings_file.cross_sections = "./mgxs.xml"
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.63, -0.63, -1, 0.63, 0.63, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC tallies.xml file
###############################################################################

# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.type = 'regular'
mesh.dimension = [100, 100, 1]
mesh.lower_left = [-0.63, -0.63, -1.e50]
mesh.upper_right = [0.63, 0.63, 1.e50]

# Instantiate some tally Filters
energy_filter = openmc.Filter(type='energy',
                              bins=[1E-11, 0.0635E-6, 10.0E-6, 1.0E-4, 1.0E-3,
                                    0.5, 1.0, 20.0])
mesh_filter = openmc.Filter()
mesh_filter.mesh = mesh

# Instantiate the Tally
tally = openmc.Tally(tally_id=1, name='tally 1')
tally.filters = [energy_filter, mesh_filter]
tally.scores = ['flux', 'fission', 'nu-fission']

# Instantiate a Tallies collection, register all Tallies, and export to XML
tallies_file = openmc.Tallies([tally])
tallies_file.export_to_xml()
