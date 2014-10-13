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
moderator.setDensity('g/cc', 1.0)
moderator.addNuclide(h1, 2.)
moderator.addNuclide(o16, 1.)
moderator.addSAlphaBeta('HH2O', '71t')

fuel = openmc.Material(material_id=40, name='fuel')
fuel.setDensity('g/cc', 4.5)
fuel.addNuclide(u235, 1.)

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.setDefaultXS('71c')
materials_file.addMaterials([moderator, fuel])
materials_file.exportToXML()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate ZCylinder surfaces
surf1 = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=7, name='surf 1')
surf2 = openmc.ZCylinder(surface_id=2, x0=0, y0=0, R=9, name='surf 2')
surf3 = openmc.ZCylinder(surface_id=3, x0=0, y0=0, R=11, name='surf 3')
surf3.setBoundaryType('vacuum')

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='cell 1')
cell2 = openmc.Cell(cell_id=100, name='cell 2')
cell3 = openmc.Cell(cell_id=101, name='cell 3')
cell4 = openmc.Cell(cell_id=2, name='cell 4')

# Register Surfaces with Cells
cell1.addSurface(surface=surf2, halfspace=-1)
cell2.addSurface(surface=surf1, halfspace=-1)
cell3.addSurface(surface=surf1, halfspace=+1)
cell4.addSurface(surface=surf2, halfspace=+1)
cell4.addSurface(surface=surf3, halfspace=-1)

# Register Materials with Cells
cell2.setFill(fuel)
cell3.setFill(moderator)
cell4.setFill(moderator)

# Instantiate Universes
universe1 = openmc.Universe(universe_id=37)
root = openmc.Universe(universe_id=0, name='root universe')
cell1.setFill(universe1)

# Register Cells with Universes
universe1.addCells([cell2, cell3])
root.addCells([cell1, cell4])

# Instantiate a Geometry and register the root Universe
geometry = openmc.Geometry()
geometry.setRootUniverse(root)

# Instantiate a GeometryFile, register Geometry, and export to XML
geometry_file = openmc.GeometryFile()
geometry_file.setGeometry(geometry)
geometry_file.exportToXML()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Instantiate a SettingsFile, set all runtime parameters, and export to XML
settings_file = openmc.SettingsFile()
settings_file.setBatches(batches)
settings_file.setInactive(inactive)
settings_file.setParticles(particles)
settings_file.setSourceSpace('box', [-4, -4, -4, 4, 4, 4])
settings_file.exportToXML()


###############################################################################
#                   Exporting to OpenMC tallies.xml File
###############################################################################

# Instantiate some tally Filters
cell_filter = openmc.Filter(type='cell', bins=100)
energy_filter = openmc.Filter(type='energy', bins=[0., 20.])
energyout_filter = openmc.Filter(type='energyout', bins=[0., 20.])

# Instantiate the first Tally
first_tally = openmc.Tally(tally_id=1, label='first tally')
first_tally.addFilter(cell_filter)
scores = ['total', 'scatter', 'nu-scatter', \
          'absorption', 'fission', 'nu-fission']
for score in scores:
  first_tally.addScore(score)

# Instantiate the second Tally
second_tally = openmc.Tally(tally_id=2, label='second tally')
second_tally.addFilter(cell_filter)
second_tally.addFilter(energy_filter)
scores = ['total', 'scatter', 'nu-scatter', \
          'absorption', 'fission', 'nu-fission']
for score in scores:
  second_tally.addScore(score)

# Instantiate the third Tally
third_tally = openmc.Tally(tally_id=3, label='third tally')
third_tally.addFilter(cell_filter)
third_tally.addFilter(energy_filter)
third_tally.addFilter(energyout_filter)
scores = ['scatter', 'nu-scatter', 'nu-fission']
for score in scores:
  third_tally.addScore(score)

# Instantiate a TalliesFile, register all Tallies, and export to XML
tallies_file = openmc.TalliesFile()
tallies_file.addTally(first_tally)
tallies_file.addTally(second_tally)
tallies_file.addTally(third_tally)
tallies_file.exportToXML()
