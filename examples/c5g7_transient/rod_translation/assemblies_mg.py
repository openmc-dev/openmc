import openmc
from pin_cells_mg import materials, surfaces, universes, cells, mgxs_lib_file

###############################################################################
#                   make abbreviates for pin cell universes
###############################################################################

###############################################################################
#                 Instantiate square assembly lattices
###############################################################################

lattices = {}

for bank in [1,4]:
    uo = universes['UO2 Bank {}'.format(bank)]
    cc = universes['Control Rod Base Bank {}'.format(bank)]
    fc = universes['Fission Chamber Bank {}'.format(bank)]
    name = 'UO2 Lattice {}'.format(bank)
    lattices[name] = openmc.RectLattice(name=name)
    lat = lattices[name]
    lat.lower_left = [-10.71, -10.71]
    lat.pitch      = [1.26, 1.26]
    lat.universes  = [[uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, uo, uo, uo],
                      [uo, uo, uo, cc, uo, uo, uo, uo, uo, uo, uo, uo, uo, cc, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, cc, uo, uo, cc, uo, uo, fc, uo, uo, cc, uo, uo, cc, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, cc, uo, uo, uo, uo, uo, uo, uo, uo, uo, cc, uo, uo, uo],
                      [uo, uo, uo, uo, uo, cc, uo, uo, cc, uo, uo, cc, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo],
                      [uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo, uo]]
    name = 'UO2 Assembly {}'.format(bank)
    universes[name]  = openmc.Universe(name=name)
    cells[name]      = openmc.Cell(name=name)
    cells[name].fill = lat
    universes[name].add_cell(cells[name])

for bank in [2,3]:
    m4 = universes['MOX 4.3% Bank {}'.format(bank)]
    m7 = universes['MOX 7.0% Bank {}'.format(bank)]
    m8 = universes['MOX 8.7% Bank {}'.format(bank)]
    cc = universes['Control Rod Base Bank {}'.format(bank)]
    fc = universes['Fission Chamber Bank {}'.format(bank)]
    name = 'MOX Lattice {}'.format(bank)
    lattices[name] = openmc.RectLattice(name=name)
    lat = lattices[name]
    lat.lower_left = [-10.71, -10.71]
    lat.pitch      = [1.26, 1.26]
    lat.universes  = [[m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4],
                      [m4, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m4],
                      [m4, m7, m7, m7, m7, cc, m7, m7, cc, m7, m7, cc, m7, m7, m7, m7, m4],
                      [m4, m7, m7, cc, m7, m8, m8, m8, m8, m8, m8, m8, m7, cc, m7, m7, m4],
                      [m4, m7, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m7, m4],
                      [m4, m7, cc, m8, m8, cc, m8, m8, cc, m8, m8, cc, m8, m8, cc, m7, m4],
                      [m4, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m4],
                      [m4, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m4],
                      [m4, m7, cc, m8, m8, cc, m8, m8, fc, m8, m8, cc, m8, m8, cc, m7, m4],
                      [m4, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m4],
                      [m4, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m4],
                      [m4, m7, cc, m8, m8, cc, m8, m8, cc, m8, m8, cc, m8, m8, cc, m7, m4],
                      [m4, m7, m7, m7, m8, m8, m8, m8, m8, m8, m8, m8, m8, m7, m7, m7, m4],
                      [m4, m7, m7, cc, m7, m8, m8, m8, m8, m8, m8, m8, m7, cc, m7, m7, m4],
                      [m4, m7, m7, m7, m7, cc, m7, m7, cc, m7, m7, cc, m7, m7, m7, m7, m4],
                      [m4, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m7, m4],
                      [m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4, m4]]
    name = 'MOX Assembly {}'.format(bank)
    universes[name]  = openmc.Universe(name=name)
    cells[name]      = openmc.Cell(name=name)
    cells[name].fill = lat
    universes[name].add_cell(cells[name])


fc = universes['Fission Chamber Bank 0']
cr = universes['Control Rod Reflector']
rf = universes['Moderator']

name = 'Control Rod Lattice'
lattices[name] = openmc.RectLattice(name=name)
lat = lattices[name]
lat.lower_left = [-10.71, -10.71]
lat.pitch      = [1.26, 1.26]
lat.universes  = [[rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, rf, rf, rf],
                  [rf, rf, rf, cr, rf, rf, rf, rf, rf, rf, rf, rf, rf, cr, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, cr, rf, rf, cr, rf, rf, fc, rf, rf, cr, rf, rf, cr, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, cr, rf, rf, rf, rf, rf, rf, rf, rf, rf, cr, rf, rf, rf],
                  [rf, rf, rf, rf, rf, cr, rf, rf, cr, rf, rf, cr, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf],
                  [rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf, rf]]
name = 'Control Rod Assembly'
universes[name]  = openmc.Universe(name=name)
cells[name]      = openmc.Cell(name=name)
cells[name].fill = lat
universes[name].add_cell(cells[name])
