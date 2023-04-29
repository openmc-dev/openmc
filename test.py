import openmc
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import openmc


# MATERIALS

breeder_material = openmc.Material()   # Pb84.2Li15.8
breeder_material.add_element('Pb', 84.2, percent_type='ao')
breeder_material.add_element('Li', 15.8, percent_type='ao', enrichment=7.0, enrichment_target='Li6', enrichment_type='ao')   # natural enrichment = 7% Li6
breeder_material.set_density('atom/b-cm', 3.2720171e-2)   # around 11 g/cm3

copper_material = openmc.Material()
copper_material.set_density('g/cm3', 8.5)
copper_material.add_element('Cu', 1.0)

eurofer_material = openmc.Material()
eurofer_material.set_density('g/cm3', 7.75)
eurofer_material.add_element('Fe', 89.067, percent_type='wo')

mats = openmc.Materials([breeder_material, eurofer_material, copper_material])

# surfaces
central_sol_surface = openmc.ZCylinder(r=100)
central_shield_outer_surface = openmc.ZCylinder(r=110)
vessel_inner_surface = openmc.Sphere(r=500)
first_wall_outer_surface = openmc.Sphere(r=510)
breeder_blanket_outer_surface = openmc.Sphere(r=610, boundary_type='vacuum')

# cells
central_sol_region = -central_sol_surface & -breeder_blanket_outer_surface
central_sol_cell = openmc.Cell(region=central_sol_region)
central_sol_cell.fill = copper_material

central_shield_region = +central_sol_surface & -central_shield_outer_surface & -breeder_blanket_outer_surface
central_shield_cell = openmc.Cell(region=central_shield_region)
central_shield_cell.fill = eurofer_material

inner_vessel_region = -vessel_inner_surface & +central_shield_outer_surface
inner_vessel_cell = openmc.Cell(region=inner_vessel_region)
# no material set as default is vacuum

first_wall_region = -first_wall_outer_surface & +vessel_inner_surface
first_wall_cell = openmc.Cell(region=first_wall_region)
first_wall_cell.fill = eurofer_material

breeder_blanket_region = +first_wall_outer_surface & -breeder_blanket_outer_surface & +central_shield_outer_surface
breeder_blanket_cell = openmc.Cell(region=breeder_blanket_region)
breeder_blanket_cell.fill = breeder_material

my_geometry = openmc.Geometry([central_sol_cell, central_shield_cell, inner_vessel_cell, first_wall_cell, breeder_blanket_cell])


# SIMULATION SETTINGS

sett = openmc.Settings()
sett.batches = 10
sett.particles = 6000
sett.run_mode = 'fixed source'

my_source = openmc.Source()
radius = openmc.stats.Discrete([200], [1])
z_values = openmc.stats.Discrete([0], [1])
angle = openmc.stats.Uniform(a=0., b=2* 3.14159265359)
my_source.space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0))
my_source.angle = openmc.stats.Isotropic()

sett.source = my_source

mesh = openmc.RegularMesh().from_domain(
    my_geometry,
    dimension=[15, 20, 25]
)

mesh_filter = openmc.MeshFilter(mesh)

mesh_tally_1 = openmc.Tally(name='tbr_on_mesh')
mesh_tally_1.filters = [mesh_filter]
mesh_tally_1.scores = ['flux']  # where X is a wildcard
tallies = openmc.Tallies([mesh_tally_1])

model = openmc.Model(my_geometry, mats, sett, tallies)
sp_filename = model.run()

sp = openmc.StatePoint(sp_filename)

tally = sp.get_tally(name='tbr_on_mesh')
tally_mean = tally.get_values()
mesh = tally.find_filter(filter_type=openmc.MeshFilter).mesh

for basis in ['xz', 'yz', 'xy']:
    slice_data = mesh.get_data_slice(
        dataset=tally_mean,
        basis = basis,
        slice_index=10  # could be changed to be middle of mesh
    )

    plt.imshow(
        slice_data,
        origin='lower',
        norm=LogNorm(),
        # extent=
    )
    
    img = my_geometry.root_universe.get_plot_image()

    rgb = (img * 256).astype(int)
    image_value = (rgb[..., 0] << 16) + \
        (rgb[..., 1] << 8) + (rgb[..., 2])
    import numpy as np
    plt.contour(
        image_value,
        origin="upper",
        colors="k",
        linestyles="solid",
        linewidths=1,
        levels=np.unique(image_value),
        # extent=(x_min, x_max, y_min, y_max),
    )
    
    plt.show()