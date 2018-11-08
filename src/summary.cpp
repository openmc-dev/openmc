#include "openmc/cell.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/surface.h"
#include "openmc/settings.h"

namespace openmc {

extern "C" void
write_geometry(hid_t file_id) {

  auto geom_group = create_group(file_id, "geometry");

#ifdef DAGMC
  if (settings::dagmc) {
    write_attribute(geom_group, "dagmc", 1);
    return;
  }
#endif

  write_attribute(geom_group, "n_cells", model::cells.size());
  write_attribute(geom_group, "n_surfaces", surfaces.size());
  write_attribute(geom_group, "n_universes", model::universes.size());
  write_attribute(geom_group, "n_lattices", model::lattices.size());

  auto cells_group = create_group(geom_group, "cells");
  for (Cell* c : model::cells) c->to_hdf5(cells_group);
  close_group(cells_group);

  auto surfaces_group = create_group(geom_group, "surfaces");
  for (Surface* surf : surfaces) surf->to_hdf5(surfaces_group);
  close_group(surfaces_group);

  auto universes_group = create_group(geom_group, "universes");
  for (Universe* u : model::universes) u->to_hdf5(universes_group);
  close_group(universes_group);

  auto lattices_group = create_group(geom_group, "lattices");
  for (Lattice* lat : model::lattices) lat->to_hdf5(lattices_group);
  close_group(lattices_group);

  close_group(geom_group);
}

} // namespace openmc
