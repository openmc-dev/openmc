import openmc
from assemblies_mg import materials, surfaces, universes, cells, lattices, mgxs_lib_file

u1 = universes['UO2 Assembly 1']
u4 = universes['UO2 Assembly 4']
m2 = universes['MOX Assembly 2']
m3 = universes['MOX Assembly 3']
r = universes['Moderator']

# Create core lattices
lattices['Core Lattice'] = openmc.RectLattice(name='Core Lattice')
lat = lattices['Core Lattice']
lat.lower_left = [-32.13, -32.13]
lat.pitch      = [ 21.42,  21.42]
lat.universes  = [[u1, m2, r],
                  [m3, u4, r],
                  [r , r , r]]

universes['Core Lattice']    = openmc.Universe(name='Core Lattice')
cells['Core Lattice']        = openmc.Cell(name='Core Lattice')
cells['Core Lattice'].fill   = lat
cells['Core Lattice'].region = +surfaces['Core x-min'] & -surfaces['Core x-max'] & \
                               +surfaces['Core y-min'] & -surfaces['Core y-max'] & \
                               +surfaces['Core z-min'] & -surfaces['Core z-max']

surfaces['Core z-min'].boundary_type = 'reflective'
surfaces['Core z-max'].boundary_type = 'reflective'

# Create core universes
universes['Core'] = openmc.Universe(universe_id=0, name='Core')
universes['Core'].add_cell(cells['Core Lattice'])

geometry = openmc.Geometry(universes['Core'])
