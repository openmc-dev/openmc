import openmc


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 500
inactive = 10
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate a Nuclides
u235 = openmc.Nuclide('U-235')

# Instantiate a Material and register the Nuclide
fuel = openmc.Material(material_id=1, name='fuel')
fuel.setDensity('g/cc', 4.5)
fuel.addNuclide(u235, 1.)

# Instantiate a MaterialsFile, register Material, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.setDefaultXS('71c')
materials_file.addMaterial(fuel)
materials_file.exportToXML()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate Surfaces
surf1 = openmc.XPlane(surface_id=1, x0=-1, name='surf 1')
surf2 = openmc.XPlane(surface_id=2, x0=+1, name='surf 2')
surf3 = openmc.YPlane(surface_id=3, y0=-1, name='surf 3')
surf4 = openmc.YPlane(surface_id=4, y0=+1, name='surf 4')
surf5 = openmc.ZPlane(surface_id=5, z0=-1, name='surf 5')
surf6 = openmc.ZPlane(surface_id=6, z0=+1, name='surf 6')

surf1.setBoundaryType('vacuum')
surf2.setBoundaryType('vacuum')
surf3.setBoundaryType('reflective')
surf4.setBoundaryType('reflective')
surf5.setBoundaryType('reflective')
surf6.setBoundaryType('reflective')

# Instantiate Cell
cell = openmc.Cell(cell_id=1, name='cell 1')

# Register Surfaces with Cell
cell.addSurface(surface=surf1, halfspace=+1)
cell.addSurface(surface=surf2, halfspace=-1)
cell.addSurface(surface=surf3, halfspace=+1)
cell.addSurface(surface=surf4, halfspace=-1)
cell.addSurface(surface=surf5, halfspace=+1)
cell.addSurface(surface=surf6, halfspace=-1)

# Register Material with Cell
cell.setFill(fuel)

# Instantiate Universes
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cell with Universe
root.addCell(cell)

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
