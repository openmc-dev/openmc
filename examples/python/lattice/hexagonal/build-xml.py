import openmc
from openmc.source import Source
from openmc.stats import Box

###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 20
inactive = 10
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')
fe56 = openmc.Nuclide('Fe-56')

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(material_id=1, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide(u235, 1.)

moderator = openmc.Material(material_id=2, name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_nuclide(h1, 2.)
moderator.add_nuclide(o16, 1.)
moderator.add_s_alpha_beta('HH2O', '71t')

iron = openmc.Material(material_id=3, name='iron')
iron.set_density('g/cc', 7.9)
iron.add_nuclide(fe56, 1.)

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.default_xs = '71c'
materials_file.add_materials([moderator, fuel, iron])
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate Surfaces
left = openmc.XPlane(surface_id=1, x0=-3, name='left')
right = openmc.XPlane(surface_id=2, x0=3, name='right')
bottom = openmc.YPlane(surface_id=3, y0=-4, name='bottom')
top = openmc.YPlane(surface_id=4, y0=4, name='top')
fuel_surf = openmc.ZCylinder(surface_id=5, x0=0, y0=0, R=0.4)

left.boundary_type = 'vacuum'
right.boundary_type = 'vacuum'
top.boundary_type = 'vacuum'
bottom.boundary_type = 'vacuum'

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')
cell2 = openmc.Cell(cell_id=101, name='cell 2')
cell3 = openmc.Cell(cell_id=102, name='cell 3')
cell4 = openmc.Cell(cell_id=500, name='cell 4')
cell5 = openmc.Cell(cell_id=600, name='cell 5')
cell6 = openmc.Cell(cell_id=601, name='cell 6')

# Use surface half-spaces to define regions
cell1.region = +left & -right & +bottom & -top
cell2.region = -fuel_surf
cell3.region = +fuel_surf
cell5.region = -fuel_surf
cell6.region = +fuel_surf

# Register Materials with Cells
cell2.fill = fuel
cell3.fill = moderator
cell4.fill = moderator
cell5.fill = iron
cell6.fill = moderator

# Instantiate Universe
univ1 = openmc.Universe(universe_id=1)
univ2 = openmc.Universe(universe_id=3)
univ3 = openmc.Universe(universe_id=4)
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
univ1.add_cells([cell2, cell3])
univ2.add_cells([cell4])
univ3.add_cells([cell5, cell6])
root.add_cell(cell1)

# Instantiate a Lattice
lattice = openmc.HexLattice(lattice_id=5)
lattice.center = [0., 0., 0.]
lattice.pitch = [1., 2.]
lattice.universes = \
    [ [ [univ2] + [univ3]*11, [univ2] + [univ3]*5, [univ3] ],
      [ [univ2] + [univ1]*11, [univ2] + [univ1]*5, [univ1] ],
      [ [univ2] + [univ3]*11, [univ2] + [univ3]*5, [univ3] ] ]
lattice.outer = univ2

# Fill Cell with the Lattice
cell1.fill = lattice

# Instantiate a Geometry and register the root Universe
geometry = openmc.Geometry()
geometry.root_universe = root

# Instantiate a GeometryFile, register Geometry, and export to XML
geometry_file = openmc.GeometryFile()
geometry_file.geometry = geometry
geometry_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Instantiate a SettingsFile, set all runtime parameters, and export to XML
settings_file = openmc.SettingsFile()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.source = Source(space=Box(
    [-1, -1, -1], [1, 1, 1]))
settings_file.keff_trigger = {'type' : 'std_dev', 'threshold' : 5E-4}
settings_file.trigger_active = True
settings_file.trigger_max_batches = 100
settings_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC plots.xml File
###############################################################################

plot_xy = openmc.Plot(plot_id=1)
plot_xy.filename = 'plot_xy'
plot_xy.origin = [0, 0, 0]
plot_xy.width = [6, 6]
plot_xy.pixels = [400, 400]
plot_xy.color = 'mat'

plot_yz = openmc.Plot(plot_id=2)
plot_yz.filename = 'plot_yz'
plot_yz.basis = 'yz'
plot_yz.origin = [0, 0, 0]
plot_yz.width = [8, 8]
plot_yz.pixels = [400, 400]
plot_yz.color = 'mat'

# Instantiate a PlotsFile, add Plot, and export to XML
plot_file = openmc.PlotsFile()
plot_file.add_plot(plot_xy)
plot_file.add_plot(plot_yz)
plot_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC tallies.xml File
###############################################################################

# Instantiate a distribcell Tally
tally = openmc.Tally(tally_id=1)
tally.filters = [openmc.Filter(type='distribcell', bins=[cell2.id])]
tally.scores = ['total']

# Instantiate a TalliesFile, register Tally/Mesh, and export to XML
tallies_file = openmc.TalliesFile()
tallies_file.add_tally(tally)
tallies_file.export_to_xml()
