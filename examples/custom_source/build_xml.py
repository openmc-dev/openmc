import openmc

# Create a single material
iron = openmc.Material()
iron.set_density('g/cm3', 5.0)
iron.add_element('Fe', 1.0)
mats = openmc.Materials([iron])
mats.export_to_xml()

# Create a 5 cm x 5 cm box filled with iron
box = openmc.model.rectangular_prism(10.0, 10.0, boundary_type='vacuum')
cell = openmc.Cell(fill=iron, region=box)
geometry = openmc.Geometry([cell])
geometry.export_to_xml()

# Tell OpenMC we're going to use our custom source
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.batches = 10
settings.particles = 1000
source = openmc.Source()
source.library = 'build/libsource.so'
settings.source = source
settings.export_to_xml()

# Finally, define a mesh tally so that we can see the resulting flux
mesh = openmc.RegularMesh()
mesh.lower_left = (-5.0, -5.0)
mesh.upper_right = (5.0, 5.0)
mesh.dimension = (50, 50)

tally = openmc.Tally()
tally.filters = [openmc.MeshFilter(mesh)]
tally.scores = ['flux']
tallies = openmc.Tallies([tally])
tallies.export_to_xml()
