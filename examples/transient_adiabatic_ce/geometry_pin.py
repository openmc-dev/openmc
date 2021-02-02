import openmc

###############################################################################
# Create materials for the problem

uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

helium = openmc.Material(name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(name='Moderator')
borated_water.set_density('g/cm3', 0.740582) 
borated_water.add_element('B',   2.7800E-5, 'ao')
borated_water.add_element('H', 2*3.3500E-2, 'ao')
borated_water.add_element('O',   3.3500E-2, 'ao')
borated_water.add_s_alpha_beta('c_H_in_H2O')

# Collect the materials together
materials = openmc.Materials([uo2, helium, zircaloy, borated_water])

###############################################################################
# Define problem geometry

# Create cylindrical surfaces
fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(r=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(r=0.45720, name='Clad OR')

# Create a region represented as the inside of a rectangular prism
pitch = 1.25984
box = openmc.rectangular_prism(pitch, pitch, boundary_type='reflective')
bottom      = openmc.ZPlane(z0=-182.88)
top         = openmc.ZPlane(z0= 182.88)
reflect_bot = openmc.ZPlane(z0=-200, boundary_type='vacuum')
reflect_top = openmc.ZPlane(z0= 200, boundary_type='vacuum')

fuel             = openmc.Cell(fill=uo2,           region=-fuel_or & +bottom & -top)
gap              = openmc.Cell(fill=helium,        region=+fuel_or & -clad_ir & +bottom & -top)
clad             = openmc.Cell(fill=zircaloy,      region=+clad_ir & -clad_or & +bottom & -top)
water            = openmc.Cell(fill=borated_water, region=+clad_or & box & +bottom & -top)
bottom_reflector = openmc.Cell(fill=borated_water, region=-bottom & +reflect_bot & box)
top_reflector    = openmc.Cell(fill=borated_water, region=+top & -reflect_top & box)

# Create a geometry
geometry = openmc.Geometry([fuel, gap, clad, water, bottom_reflector, top_reflector])
