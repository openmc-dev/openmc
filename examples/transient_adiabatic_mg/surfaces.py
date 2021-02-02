import openmc

###############################################################################
#                       Create the OpenMC Surfaces
###############################################################################

surfaces = {}

# Instantiate ZCylinder surfaces
surfaces['Fuel OR']            = openmc.ZCylinder(x0=0, y0=0, R=0.4095, name='Fuel OR')
surfaces['Fuel Inner Clad IR'] = openmc.ZCylinder(x0=0, y0=0, R=0.4180, name='Fuel Inner Clad IR')
surfaces['Fuel Inner Clad OR'] = openmc.ZCylinder(x0=0, y0=0, R=0.4750, name='Fuel Inner Clad OR')
surfaces['Fuel Outer Clad IR'] = openmc.ZCylinder(x0=0, y0=0, R=0.4800, name='Fuel Outer Clad IR')
surfaces['Fuel Outer Clad OR'] = openmc.ZCylinder(x0=0, y0=0, R=0.5400, name='Fuel Outer Clad OR')
surfaces['Guide Tube IR']      = openmc.ZCylinder(x0=0, y0=0, R=0.3400, name='Guide Tube IR')
surfaces['Guide Tube OR']      = openmc.ZCylinder(x0=0, y0=0, R=0.5400, name='Guide Tube IR')

# Instantiate the axial surfaces
surfaces['Axial Midplane']        = openmc.ZPlane(z0=0.0   , name='Axial Midplane')
surfaces['Core x-min']            = openmc.XPlane(x0=-32.13, name='Core x-min')
surfaces['Core x-max']            = openmc.XPlane(x0= 32.13, name='Core x-max')
surfaces['Core y-min']            = openmc.YPlane(y0=-32.13, name='Core y-min')
surfaces['Core y-max']            = openmc.YPlane(y0= 32.13, name='Core y-max')
surfaces['Core z-min']            = openmc.ZPlane(z0=-64.26, name='Core z-min')
surfaces['Core z-max']            = openmc.ZPlane(z0= 64.26, name='Core z-max')
surfaces['Lower Reflector z-min'] = openmc.ZPlane(z0=-85.68, name='Lower Reflector z-min')
surfaces['Upper Reflector z-max'] = openmc.ZPlane(z0= 85.68, name='Upper Reflector z-max')


# Set the boundary conditions
surfaces['Core x-min'].boundary_type = 'reflective'
surfaces['Core x-max'].boundary_type = 'vacuum'
surfaces['Core y-min'].boundary_type = 'vacuum'
surfaces['Core y-max'].boundary_type = 'reflective'
surfaces['Lower Reflector z-min'].boundary_type = 'reflective'
surfaces['Upper Reflector z-max'].boundary_type = 'reflective'
