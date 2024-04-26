import openmc
import matplotlib as plt

universe_1 = openmc.Universe(universe_id=4, name='root universe')
universe_2 = openmc.Universe(universe_id=5, name='new universe')
latt = openmc.RectLattice()
latt.pitch = (4,4)
lattice_fill = [[universe_1, universe_2]]
latt.universes = lattice_fill

parent_cell = openmc.Cell(fill=latt)
parent_cell.plot()
plt.show()
