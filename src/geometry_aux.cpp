#include "openmc/geometry_aux.h"

#include <algorithm>  // for std::max
#include <sstream>
#include <unordered_set>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/lattice.h"
#include "openmc/material.h"


namespace openmc {

//==============================================================================

void
adjust_indices()
{
  // Adjust material/fill idices.
  for (Cell* c : global_cells) {
    if (c->fill != C_NONE) {
      int32_t id = c->fill;
      auto search_univ = universe_map.find(id);
      auto search_lat = lattice_map.find(id);
      if (search_univ != universe_map.end()) {
        c->type = FILL_UNIVERSE;
        c->fill = search_univ->second;
      } else if (search_lat != lattice_map.end()) {
        c->type = FILL_LATTICE;
        c->fill = search_lat->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Specified fill " << id << " on cell " << c->id
                << " is neither a universe nor a lattice.";
        fatal_error(err_msg);
      }
    } else {
      c->type = FILL_MATERIAL;
      for (auto it = c->material.begin(); it != c->material.end(); it++) {
        int32_t mid = *it;
        if (mid != MATERIAL_VOID) {
          auto search = material_map.find(mid);
          if (search != material_map.end()) {
            *it = search->second;
          } else {
            std::stringstream err_msg;
            err_msg << "Could not find material " << mid
                    << " specified on cell " << c->id;
            fatal_error(err_msg);
          }
        }
      }
    }
  }

  // Change cell.universe values from IDs to indices.
  for (Cell* c : global_cells) {
    auto search = universe_map.find(c->universe);
    if (search != universe_map.end()) {
      //TODO: Remove this off-by-one indexing.
      c->universe = search->second + 1;
    } else {
      std::stringstream err_msg;
      err_msg << "Could not find universe " << c->universe
              << " specified on cell " << c->id;
      fatal_error(err_msg);
    }
  }

  // Change all lattice universe values from IDs to indices.
  for (Lattice* l : lattices_c) {
    l->adjust_indices();
  }
}

//==============================================================================

int32_t
find_root_universe()
{
  // Find all the universes listed as a cell fill.
  std::unordered_set<int32_t> fill_univ_ids;
  for (Cell* c : global_cells) {
    fill_univ_ids.insert(c->fill);
  }

  // Find all the universes contained in a lattice.
  for (Lattice* lat : lattices_c) {
    for (auto it = lat->begin(); it != lat->end(); ++it) {
      fill_univ_ids.insert(*it);
    }
    if (lat->outer != NO_OUTER_UNIVERSE) {
      fill_univ_ids.insert(lat->outer);
    }
  }

  // Figure out which universe is not in the set.  This is the root universe.
  bool root_found {false};
  int32_t root_univ;
  for (int32_t i = 0; i < global_universes.size(); i++) {
    auto search = fill_univ_ids.find(global_universes[i]->id);
    if (search == fill_univ_ids.end()) {
      if (root_found) {
        fatal_error("Two or more universes are not used as fill universes, so "
                    "it is not possible to distinguish which one is the root "
                    "universe.");
      } else {
        root_found = true;
        root_univ = i;
      }
    }
  }
  if (!root_found) fatal_error("Could not find a root universe.  Make sure "
       "there are no circular dependencies in the geometry.");

  return root_univ;
}

//==============================================================================

void
allocate_offset_tables(int n_maps)
{
  for (Cell* c : global_cells) {
    if (c->type != FILL_MATERIAL) {
      c->offset.resize(n_maps, C_NONE);
    }
  }

  for (Lattice* lat : lattices_c) {
    lat->allocate_offset_table(n_maps);
  }
}

//==============================================================================

void
count_cell_instances(int32_t univ_indx)
{
  for (int32_t cell_indx : global_universes[univ_indx]->cells) {
    Cell& c = *global_cells[cell_indx];
    ++c.n_instances;

    if (c.type == FILL_UNIVERSE) {
      // This cell contains another universe.  Recurse into that universe.
      count_cell_instances(c.fill);

    } else if (c.type == FILL_LATTICE) {
      // This cell contains a lattice.  Recurse into the lattice universes.
      Lattice& lat = *lattices_c[c.fill];
      for (auto it = lat.begin(); it != lat.end(); ++it) {
        count_cell_instances(*it);
      }
    }
  }
}

//==============================================================================

int
count_universe_instances(int32_t search_univ, int32_t target_univ_id)
{
  //  If this is the target, it can't contain itself.
  if (global_universes[search_univ]->id == target_univ_id) {
    return 1;
  }

  int count {0};
  for (int32_t cell_indx : global_universes[search_univ]->cells) {
    Cell& c = *global_cells[cell_indx];

    if (c.type == FILL_UNIVERSE) {
      int32_t next_univ = c.fill;
      count += count_universe_instances(next_univ, target_univ_id);

    } else if (c.type == FILL_LATTICE) {
      Lattice& lat = *lattices_c[c.fill];
      for (auto it = lat.begin(); it != lat.end(); ++it) {
        int32_t next_univ = *it;
        count += count_universe_instances(next_univ, target_univ_id);
      }
    }
  }

  return count;
}

//==============================================================================

void
fill_offset_tables(int32_t target_univ_id, int map)
{
  for (Universe* univ : global_universes) {
    int32_t offset {0};  // TODO: is this a bug?  It matches F90 implementation.
    for (int32_t cell_indx : univ->cells) {
      Cell& c = *global_cells[cell_indx];

      if (c.type == FILL_UNIVERSE) {
        c.offset[map] = offset;
        int32_t search_univ = c.fill;
        offset += count_universe_instances(search_univ, target_univ_id);

      } else if (c.type == FILL_LATTICE) {
        Lattice& lat = *lattices_c[c.fill];
        offset = lat.fill_offset_table(offset, target_univ_id, map);
      }
    }
  }
}

//==============================================================================

std::string
distribcell_path_inner(int32_t target_cell, int32_t map, int32_t target_offset,
                       const Universe& search_univ, int32_t offset)
{
  std::stringstream path;

  path << "u" << search_univ.id << "->";

  // Check to see if this universe directly contains the target cell.  If so,
  // write to the path and return.
  for (int32_t cell_indx : search_univ.cells) {
    if ((cell_indx == target_cell) && (offset == target_offset)) {
      Cell& c = *global_cells[cell_indx];
      path << "c" << c.id;
      return path.str();
    }
  }

  // The target must be further down the geometry tree and contained in a fill
  // cell or lattice cell in this universe.  Find which cell contains the
  // target.
  std::vector<std::int32_t>::const_reverse_iterator cell_it
       {search_univ.cells.crbegin()};
  for (; cell_it != search_univ.cells.crend(); ++cell_it) {
    Cell& c = *global_cells[*cell_it];

    // Material cells don't contain other cells so ignore them.
    if (c.type != FILL_MATERIAL) {
      int32_t temp_offset;
      if (c.type == FILL_UNIVERSE) {
        temp_offset = offset + c.offset[map];
      } else {
        Lattice& lat = *lattices_c[c.fill];
        int32_t indx = lat.universes.size()*map + lat.begin().indx;
        temp_offset = offset + lat.offsets[indx];
      }

      // The desired cell is the first cell that gives an offset smaller or
      // equal to the target offset.
      if (temp_offset <= target_offset) break;
    }
  }

  // Add the cell to the path string.
  Cell& c = *global_cells[*cell_it];
  path << "c" << c.id << "->";

  if (c.type == FILL_UNIVERSE) {
    // Recurse into the fill cell.
    offset += c.offset[map];
    path << distribcell_path_inner(target_cell, map, target_offset,
                                   *global_universes[c.fill], offset);
    return path.str();
  } else {
    // Recurse into the lattice cell.
    Lattice& lat = *lattices_c[c.fill];
    path << "l" << lat.id;
    for (ReverseLatticeIter it = lat.rbegin(); it != lat.rend(); ++it) {
      int32_t indx = lat.universes.size()*map + it.indx;
      int32_t temp_offset = offset + lat.offsets[indx];
      if (temp_offset <= target_offset) {
        offset = temp_offset;
        path << "(" << lat.index_to_string(it.indx) << ")->";
        path << distribcell_path_inner(target_cell, map, target_offset,
                                       *global_universes[*it], offset);
        return path.str();
      }
    }
  }
}

//==============================================================================

int
distribcell_path_len(int32_t target_cell, int32_t map, int32_t target_offset,
                     int32_t root_univ)
{
  Universe& root = *global_universes[root_univ];
  std::string path_ {distribcell_path_inner(target_cell, map, target_offset,
                                            root, 0)};
  return path_.size() + 1;
}

//==============================================================================

void
distribcell_path(int32_t target_cell, int32_t map, int32_t target_offset,
                 int32_t root_univ, char* path)
{
  Universe& root = *global_universes[root_univ];
  std::string path_ {distribcell_path_inner(target_cell, map, target_offset,
                                            root, 0)};
  path_.copy(path, path_.size());
  path[path_.size()] = '\0';
}

//==============================================================================

int
maximum_levels(int32_t univ)
{
  int levels_below {0};

  for (int32_t cell_indx : global_universes[univ]->cells) {
    Cell& c = *global_cells[cell_indx];
    if (c.type == FILL_UNIVERSE) {
      int32_t next_univ = c.fill;
      levels_below = std::max(levels_below, maximum_levels(next_univ));
    } else if (c.type == FILL_LATTICE) {
      Lattice& lat = *lattices_c[c.fill];
      for (auto it = lat.begin(); it != lat.end(); ++it) {
        int32_t next_univ = *it;
        levels_below = std::max(levels_below, maximum_levels(next_univ));
      }
    }
  }

  ++levels_below;
  return levels_below;
}

//==============================================================================

void
free_memory_geometry_c()
{
  for (Cell* c : global_cells) {delete c;}
  global_cells.clear();
  cell_map.clear();
  n_cells = 0;

  for (Universe* u : global_universes) {delete u;}
  global_universes.clear();
  universe_map.clear();

  for (Lattice* lat : lattices_c) {delete lat;}
  lattices_c.clear();
  lattice_map.clear();
}

} // namespace openmc
