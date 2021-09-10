import openmc
from assemblies_mg import materials, surfaces, universes, cells, lattices, mgxs_lib_file

u1 = universes['UO2 Assembly 1']
u4 = universes['UO2 Assembly 4']
m2 = universes['MOX Assembly 2']
m3 = universes['MOX Assembly 3']
c = universes['Control Rod Assembly']
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

lattices['Upper Reflector Lattice'] = openmc.RectLattice(name='Upper Reflector Lattice')
lat = lattices['Upper Reflector Lattice']
lat.lower_left = [-32.13, -32.13]
lat.pitch      = [ 21.42,  21.42]
lat.universes  = [[c, c, r],
                  [c, c, r],
                  [r, r, r]]
universes['Upper Reflector Lattice']    = openmc.Universe(name='Upper Reflector Lattice')
cells['Upper Reflector Lattice']        = openmc.Cell(name='Upper Reflector Lattice')
cells['Upper Reflector Lattice'].fill   = lat
cells['Upper Reflector Lattice'].region = +surfaces['Core x-min'] & -surfaces['Core x-max'] & \
                                          +surfaces['Core y-min'] & -surfaces['Core y-max'] & \
                                          +surfaces['Core z-max'] & -surfaces['Upper Reflector z-max']

cells['Lower Reflector']        = openmc.Cell(name='Lower Reflector')
cells['Lower Reflector'].fill   = materials['Moderator']
cells['Lower Reflector'].region =  +surfaces['Core x-min'] & -surfaces['Core x-max'] & \
                                   +surfaces['Core y-min'] & -surfaces['Core y-max'] & \
                                   -surfaces['Core z-min'] & +surfaces['Lower Reflector z-min']

# Create core universes
universes['Core'] = openmc.Universe(universe_id=0, name='Core')
universes['Core'].add_cells([cells['Core Lattice'],
                             cells['Upper Reflector Lattice'],
                             cells['Lower Reflector']])

geometry = openmc.Geometry(universes['Core'])
