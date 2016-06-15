import openmc


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 15
inactive = 5
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')

# Instantiate some Materials and register the appropriate Nuclides
moderator = openmc.Material(material_id=41, name='moderator')
moderator.set_density('g/cc', 1.0)
moderator.add_nuclide(h1, 2.)
moderator.add_nuclide(o16, 1.)
moderator.add_s_alpha_beta('HH2O', '71t')

fuel = openmc.Material(material_id=40, name='fuel')
fuel.set_density('g/cc', 4.5)
fuel.add_nuclide(u235, 1.)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([moderator, fuel])
materials_file.default_xs = '71c'
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate ZCylinder surfaces
surf1 = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=7, name='surf 1')
surf2 = openmc.ZCylinder(surface_id=2, x0=0, y0=0, R=9, name='surf 2')
surf3 = openmc.ZCylinder(surface_id=3, x0=0, y0=0, R=11, name='surf 3')
surf3.boundary_type = 'vacuum'

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='cell 1')
cell2 = openmc.Cell(cell_id=100, name='cell 2')
cell3 = openmc.Cell(cell_id=101, name='cell 3')
cell4 = openmc.Cell(cell_id=2, name='cell 4')

# Use surface half-spaces to define regions
cell1.region = -surf2
cell2.region = -surf1
cell3.region = +surf1
cell4.region = +surf2 & -surf3

# Register Materials with Cells
cell2.fill = fuel
cell3.fill = moderator
cell4.fill = moderator

# Instantiate Universes
universe1 = openmc.Universe(universe_id=37)
root = openmc.Universe(universe_id=0, name='root universe')
cell1.fill = universe1

# Register Cells with Universes
universe1.add_cells([cell2, cell3])
root.add_cells([cell1, cell4])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-4., -4., -4., 4., 4., 4.]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

settings_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC tallies.xml file
###############################################################################

# Instantiate some tally Filters
cell_filter = openmc.Filter(type='cell', bins=100)
energy_filter = openmc.Filter(type='energy', bins=[0., 20.])
energyout_filter = openmc.Filter(type='energyout', bins=[0., 20.])

# Instantiate the first Tally
first_tally = openmc.Tally(tally_id=1, name='first tally')
first_tally.filters = [cell_filter]
scores = ['total', 'scatter', 'nu-scatter',
          'absorption', 'fission', 'nu-fission']
first_tally.scores = scores

# Instantiate the second Tally
second_tally = openmc.Tally(tally_id=2, name='second tally')
second_tally.filters = [cell_filter, energy_filter]
second_tally.scores = scores

# Instantiate the third Tally
third_tally = openmc.Tally(tally_id=3, name='third tally')
third_tally.filters = [cell_filter, energy_filter, energyout_filter]
third_tally.scores = ['scatter', 'nu-scatter', 'nu-fission']

# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies((first_tally, second_tally, third_tally))
tallies_file.export_to_xml()
