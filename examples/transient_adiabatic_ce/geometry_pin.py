import openmc

###############################################################################
#                Create the OpenMC materials
###############################################################################

uo2 = openmc.Material(material_id=1, name='UO2')
uo2.set_density('g/cm3', 10.29769)
uo2.add_nuclide('U235', 8.6500E-4 , 'ao') 
uo2.add_nuclide('U238', 2.2250E-2 , 'ao') 
uo2.add_element('O'   , 4.62200E-2, 'ao')

water = openmc.Material(name='Moderator')
water.set_density('g/cm3', 0.74) 
water.add_element('H', 2*3.3500E-2, 'ao')
water.add_element('O',   3.3500E-2, 'ao')
water.add_element('B',   2.7800E-5, 'ao')
water.add_s_alpha_beta('c_H_in_H2O')

zirc = openmc.Material(material_id=3, name='Zr Clad')
zirc.set_density('g/cm3', 6.55)
zirc.add_element('Zr', 4.3000E-2, 'ao')

void = openmc.Material(material_id=4, name='Void')
void.set_density('g/cm3', 0.001598)
void.add_element('He', 1.e-10, 'ao')

materials = openmc.Materials([uo2, water, zirc, void])

###############################################################################
#                       Create the OpenMC Surfaces
###############################################################################

# Instantiate ZCylinder surfaces
fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(x0=0, y0=0, r=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.4572, name='Clad OR')

# Instantiate the axial surfaces
left = openmc.XPlane(x0=-0.62992, name='Pin x-min', boundary_type='reflective')
right = openmc.XPlane(x0= 0.62992, name='Pin x-max', boundary_type='reflective')
back = openmc.YPlane(y0=-0.62992, name='Pin y-min', boundary_type='reflective')
front = openmc.YPlane(y0= 0.62992, name='Pin y-max', boundary_type='reflective')
bottom = openmc.ZPlane(z0=-182.88, name='Pin z-min')
top = openmc.ZPlane(z0= 182.88, name='Pin z-max')
reflect_bot = openmc.ZPlane(z0=-200, name='Lower Reflector z-min', boundary_type='vacuum')
reflect_top = openmc.ZPlane(z0= 200, name='Upper Reflector z-max', boundary_type='vacuum')

###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

fuel = openmc.Cell(name='Fuel')
fuel.region = -fuel_or & +bottom & -top
fuel.fill = uo2

gap = openmc.Cell(name='Gap')
gap.region = +fuel_or & -clad_ir & +bottom & -top
gap.fill = void

clad = openmc.Cell(name='Clad')
clad.region = +clad_ir & -clad_or & +bottom & -top
clad.fill = zirc

moderator = openmc.Cell(name='Moderator')
moderator.region = +clad_or & -right & +left & -front & +back & +bottom & -top
moderator.fill = water

BottomReflector = openmc.Cell(name='BottomReflector')
BottomReflector.region = -bottom & +reflect_bot & -right & +left & -front & +back
BottomReflector.fill = water

TopReflector = openmc.Cell(name='TopReflector')
TopReflector.region = +top & -reflect_top & -right & +left & -front & +back
TopReflector.fill = water

root = openmc.Universe(universe_id=0, name='Root')
root.add_cells([fuel, gap, clad, moderator, BottomReflector, TopReflector])
geometry = openmc.Geometry(root)
