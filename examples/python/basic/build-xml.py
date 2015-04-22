import openmc


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 15
inactive = 5
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
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

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.default_xs = '71c'
materials_file.add_materials([moderator, fuel])
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
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

# Register Surfaces with Cells
cell1.add_surface(surface=surf2, halfspace=-1)
cell2.add_surface(surface=surf1, halfspace=-1)
cell3.add_surface(surface=surf1, halfspace=+1)
cell4.add_surface(surface=surf2, halfspace=+1)
cell4.add_surface(surface=surf3, halfspace=-1)

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
settings_file.set_source_space('box', [-4, -4, -4, 4, 4, 4])
settings_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC tallies.xml File
###############################################################################

# Instantiate some tally Filters
cell_filter = openmc.Filter(type='cell', bins=100)
energy_filter = openmc.Filter(type='energy', bins=[0., 20.])
energyout_filter = openmc.Filter(type='energyout', bins=[0., 20.])

# Instantiate the first Tally
first_tally = openmc.Tally(tally_id=1, label='first tally')
first_tally.add_filter(cell_filter)
scores = ['total', 'scatter', 'nu-scatter', \
          'absorption', 'fission', 'nu-fission']
for score in scores:
  first_tally.add_score(score)

# Instantiate the second Tally
second_tally = openmc.Tally(tally_id=2, label='second tally')
second_tally.add_filter(cell_filter)
second_tally.add_filter(energy_filter)
scores = ['total', 'scatter', 'nu-scatter', \
          'absorption', 'fission', 'nu-fission']
for score in scores:
  second_tally.add_score(score)

# Instantiate the third Tally
third_tally = openmc.Tally(tally_id=3, label='third tally')
third_tally.add_filter(cell_filter)
third_tally.add_filter(energy_filter)
third_tally.add_filter(energyout_filter)
scores = ['scatter', 'nu-scatter', 'nu-fission']
for score in scores:
  third_tally.add_score(score)

# Instantiate a TalliesFile, register all Tallies, and export to XML
tallies_file = openmc.TalliesFile()
tallies_file.add_tally(first_tally)
tallies_file.add_tally(second_tally)
tallies_file.add_tally(third_tally)
tallies_file.export_to_xml()
