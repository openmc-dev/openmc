import openmc
import numpy as np
import openmc.model
import openmc.stats

fuel = openmc.Material()
fuel.add_element('U',1.0,enrichment=4.5)
fuel.add_element('O',2.0)
fuel.set_density('g/cm3',10.0)

water = openmc.Material()
water.add_element('H',2.0)
water.add_element('O',1.0)
water.set_density('g/cm3',0.74)
water.add_s_alpha_beta('c_H_in_H2O')

mats = openmc.Materials([fuel, water])

r_fuel=0.392
N_rings=10

ring_radii=[]
for i in range(N_rings + 1):
    ring_radii.append(r_fuel*np.sqrt(i/N_rings))

pitch=1.26

fuel_surf = openmc.ZCylinder(r=r_fuel)
#clad_surf = openmc.ZCylinder(r=0.46)
box= openmc.stats.Box([-pitch/2,-pitch/2,-1],[pitch/2,pitch/2,1],only_fissionable=True)

minx=openmc.XPlane(x0=-pitch/2,boundary_type='reflective')
maxx=openmc.XPlane(x0=pitch/2,boundary_type='reflective')
miny=openmc.YPlane(y0=-pitch/2,boundary_type='reflective')
maxy=openmc.YPlane(y0=pitch/2,boundary_type='reflective')

water_region= +fuel_surf & + minx & -maxx & +miny & -maxy

fuel_cell=openmc.Cell(region=-fuel_surf, fill=fuel)
water_cell=openmc.Cell(region=water_region, fill=water)

geom=openmc.Geometry([fuel_cell,water_cell])

mesh = openmc.CylindricalMesh(r_grid=ring_radii,z_grid=[-1,1],phi_grid=[0,2*np.pi])

tally_truth = openmc.Tally(name='tally_truth')
tally_truth.filters = [openmc.MeshFilter(mesh)]
tally_truth.scores = ['absorption']
tally_truth.nuclides = ['U238']

zernike_filter= openmc.ZernikeRadialFilter(order=20,r=r_fuel)

tally_zernike=openmc.Tally(name='zernike_tally')
tally_zernike.filters= [zernike_filter]
tally_zernike.scores=['absorption']
tally_zernike.nuclides=['U238']


tallies = openmc.Tallies([tally_truth, tally_zernike])


settings=openmc.Settings()
settings.batches=50
settings.inactive=10
settings.particles=5000
settings.source= openmc.Source(space=box)

model=openmc.Model(materials=mats, geometry=geom, settings=settings, tallies=tallies)

if __name__=="__main__":
    model.run()