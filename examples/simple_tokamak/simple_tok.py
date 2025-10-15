import openmc
import numpy as np

model = openmc.Model.from_model_xml(path="original_mc_model.xml")

# Adjust original MC boundaries to fit more tightly for use with random ray
reflective_face_1 = model.geometry.get_all_surfaces()[10000]
reflective_face_2 = model.geometry.get_all_surfaces()[11000]

vac_zmax = openmc.ZPlane(z0=  1000, boundary_type='vacuum')
vac_zmin = openmc.ZPlane(z0= -1000, boundary_type='vacuum')
vac_ymax = openmc.YPlane(y0=  625, boundary_type='vacuum')
vac_ymin = openmc.YPlane(y0= -625, boundary_type='vacuum')
vac_xmax = openmc.XPlane(x0= 1864.0, boundary_type='vacuum')
vac_xmin = openmc.XPlane(x0= -100.0, boundary_type='vacuum')
                
outer_region = (
    -vac_zmax & +vac_zmin &
    -vac_ymax & +vac_ymin &
    -vac_xmax & +vac_xmin &
    +reflective_face_1 & -reflective_face_2
)

outer_cell = openmc.Cell(region=outer_region, 
                            fill=model.geometry.root_universe)
outer_universe  = openmc.Universe(cells=[outer_cell])
model.geometry.root_universe = outer_universe
model.geometry.determine_paths()

# Make coarse mesh and tally for random ray void SR subdivision
mesh_cell_size_cm = 100.0 # cm
coarse_mesh = openmc.RegularMesh()
bbox = model.geometry.bounding_box
ll   = np.array(bbox.lower_left)
ur   = np.array(bbox.upper_right)
coarse_mesh.lower_left  = ll
coarse_mesh.upper_right = ur
dims = np.ceil((ur - ll) / mesh_cell_size_cm).astype(int)
coarse_mesh.dimension   = tuple(dims)
coarse_mesh_filter = openmc.MeshFilter(coarse_mesh)
coarse_mesh_tally = openmc.Tally(name="coarse_mesh_tally")
coarse_mesh_tally.filters = [coarse_mesh_filter]
coarse_mesh_tally.scores = ["flux"]
model.tallies.append(coarse_mesh_tally)


# Make fine mesh and tally for random ray void SR subdivision
mesh_cell_size_cm = 7.5 # cm
fine_mesh = openmc.RegularMesh()
bbox = model.geometry.bounding_box
ll   = np.array(bbox.lower_left)
ur   = np.array(bbox.upper_right)
fine_mesh.lower_left  = ll
fine_mesh.upper_right = ur
dims = np.ceil((ur - ll) / mesh_cell_size_cm).astype(int)
fine_mesh.dimension   = tuple(dims)
fine_mesh_filter = openmc.MeshFilter(fine_mesh)
fine_mesh_tally = openmc.Tally(name="fine_mesh_tally")
fine_mesh_tally.filters = [fine_mesh_filter]
fine_mesh_tally.scores = ["flux"]
model.tallies.append(fine_mesh_tally)

model.convert_to_multigroup(
    method = "stochastic_slab",
    nparticles = 10000,
    groups='CASMO-8',
    correction=None,
    overwrite_mgxs_library=False
)

model.convert_to_random_ray()  

bbox   = model.geometry.bounding_box
orig_src = model.settings.source[0]

rr_src = openmc.IndependentSource()
rr_src.space       = openmc.stats.Box(bbox.lower_left, bbox.upper_right)
rr_src.constraints = {'domains': [outer_cell]}

void_cells = [c for c in model.geometry.get_all_cells().values() if c.fill is None]
material_cells = [
    c for c in model.geometry.get_all_cells().values()
    if isinstance(c.fill, openmc.Material)
]

model.settings.random_ray['ray_source'] = rr_src
model.settings.random_ray["source_region_meshes"] = [(fine_mesh, material_cells), (coarse_mesh, void_cells)]
model.settings.random_ray['volume_estimator'] = 'naive'
model.settings.random_ray["distance_inactive"] = 1500.0
model.settings.random_ray["distance_active"] = 3000.0
model.settings.random_ray["sample_method"] = 'prng'
model.settings.particles = 40000
model.settings.batches   = 200
model.settings.inactive  = 100

plot = openmc.Plot()
box = model.geometry.bounding_box
plot.origin = (0.5*(box.lower_left[0]+box.upper_right[0]),
               0.5*(box.lower_left[1]+box.upper_right[1]),
               0.5*(box.lower_left[2]+box.upper_right[2]))
plot.width = [box.upper_right[0]-box.lower_left[0],
              box.upper_right[1]-box.lower_left[1],
              box.upper_right[2]-box.lower_left[2]]
plot.pixels = [300, 300, 300]
plot.type = 'voxel'
model.plots = [plot]

model.export_to_model_xml()