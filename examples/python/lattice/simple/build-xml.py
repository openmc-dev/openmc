import openmc


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

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(material_id=1, name='fuel')
fuel.setDensity('g/cc', 4.5)
fuel.addNuclide(u235, 1.)

moderator = openmc.Material(material_id=2, name='moderator')
moderator.setDensity('g/cc', 1.0)
moderator.addNuclide(h1, 2.)
moderator.addNuclide(o16, 1.)
moderator.addSAlphaBeta('HH2O', '71t')

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.setDefaultXS('71c')
materials_file.addMaterials([moderator, fuel])
materials_file.exportToXML()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate Surfaces
left = openmc.XPlane(surface_id=1, x0=-2, name='left')
right = openmc.XPlane(surface_id=2, x0=2, name='right')
bottom = openmc.YPlane(surface_id=3, y0=-2, name='bottom')
top = openmc.YPlane(surface_id=4, y0=2, name='top')
fuel1 = openmc.ZCylinder(surface_id=5, x0=0, y0=0, R=0.4)
fuel2 = openmc.ZCylinder(surface_id=6, x0=0, y0=0, R=0.3)
fuel3 = openmc.ZCylinder(surface_id=7, x0=0, y0=0, R=0.2)

left.setBoundaryType('vacuum')
right.setBoundaryType('vacuum')
top.setBoundaryType('vacuum')
bottom.setBoundaryType('vacuum')

# Instantiate Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')
cell2 = openmc.Cell(cell_id=101, name='cell 2')
cell3 = openmc.Cell(cell_id=102, name='cell 3')
cell4 = openmc.Cell(cell_id=201, name='cell 4')
cell5 = openmc.Cell(cell_id=202, name='cell 5')
cell6 = openmc.Cell(cell_id=301, name='cell 6')
cell7 = openmc.Cell(cell_id=302, name='cell 7')

# Register Surfaces with Cells
cell1.addSurface(left, halfspace=+1)
cell1.addSurface(right, halfspace=-1)
cell1.addSurface(bottom, halfspace=+1)
cell1.addSurface(top, halfspace=-1)
cell2.addSurface(fuel1, halfspace=-1)
cell3.addSurface(fuel1, halfspace=+1)
cell4.addSurface(fuel2, halfspace=-1)
cell5.addSurface(fuel2, halfspace=+1)
cell6.addSurface(fuel3, halfspace=-1)
cell7.addSurface(fuel3, halfspace=+1)

# Register Materials with Cells
cell2.setFill(fuel)
cell3.setFill(moderator)
cell4.setFill(fuel)
cell5.setFill(moderator)
cell6.setFill(fuel)
cell7.setFill(moderator)

# Instantiate Universe
univ1 = openmc.Universe(universe_id=1)
univ2 = openmc.Universe(universe_id=2)
univ3 = openmc.Universe(universe_id=3)
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
univ1.addCells([cell2, cell3])
univ2.addCells([cell4, cell5])
univ3.addCells([cell6, cell7])
root.addCell(cell1)

# Instantiate a Lattice
lattice = openmc.Lattice(lattice_id=5)
lattice.setDimension([4, 4])
lattice.setLowerLeft([-2., -2.])
lattice.setWidth([1., 1.])
lattice.setUniverses([[univ1, univ2, univ1, univ2],
                      [univ2, univ3, univ2, univ3],
                      [univ1, univ2, univ1, univ2],
                      [univ2, univ3, univ2, univ3]])

# Fill Cell with the Lattice
cell1.setFill(lattice)

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
settings_file.setSourceSpace('box', [-1, -1, -1, 1, 1, 1])
settings_file.exportToXML()


###############################################################################
#                   Exporting to OpenMC plots.xml File
###############################################################################

plot = openmc.Plot(plot_id=1)
plot.setOrigin([0, 0, 0])
plot.setWidth([4, 4])
plot.setPixels([400, 400])
plot.setColor('mat')

# Instantiate a PlotsFile, add Plot, and export to XML
plot_file = openmc.PlotsFile()
plot_file.addPlot(plot)
plot_file.exportToXML()


###############################################################################
#                   Exporting to OpenMC tallies.xml File
###############################################################################

# Instantiate a tally mesh
mesh = openmc.Mesh(mesh_id=1)
mesh.setType('rectangular')
mesh.setDimension([4, 4])
mesh.setLowerLeft([-2, -2])
mesh.setWidth([1, 1])

# Instantiate tally Filter
mesh_filter = openmc.Filter()
mesh_filter.setMesh(mesh)

# Instantiate the Tally
tally = openmc.Tally(tally_id=1)
tally.addFilter(mesh_filter)
tally.addScore('total')

# Instantiate a TalliesFile, register Tally/Mesh, and export to XML
tallies_file = openmc.TalliesFile()
tallies_file.addMesh(mesh)
tallies_file.addTally(tally)
tallies_file.exportToXML()
