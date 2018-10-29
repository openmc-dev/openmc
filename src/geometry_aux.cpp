#include "openmc/geometry_aux.h"

#include <algorithm>  // for std::max
#include <sstream>
#include <unordered_set>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_distribcell.h"


namespace openmc {

//==============================================================================

void
adjust_indices()
{
  // Adjust material/fill idices.
  for (Cell* c : cells) {
    if (c->fill_ != C_NONE) {
      int32_t id = c->fill_;
      auto search_univ = universe_map.find(id);
      auto search_lat = lattice_map.find(id);
      if (search_univ != universe_map.end()) {
        c->type_ = FILL_UNIVERSE;
        c->fill_ = search_univ->second;
      } else if (search_lat != lattice_map.end()) {
        c->type_ = FILL_LATTICE;
        c->fill_ = search_lat->second;
      } else {
        std::stringstream err_msg;
        err_msg << "Specified fill " << id << " on cell " << c->id_
                << " is neither a universe nor a lattice.";
        fatal_error(err_msg);
      }
    } else {
      c->type_ = FILL_MATERIAL;
      for (auto it = c->material_.begin(); it != c->material_.end(); it++) {
        int32_t mid = *it;
        if (mid != MATERIAL_VOID) {
          auto search = material_map.find(mid);
          if (search != material_map.end()) {
            *it = search->second;
          } else {
            std::stringstream err_msg;
            err_msg << "Could not find material " << mid
                    << " specified on cell " << c->id_;
            fatal_error(err_msg);
          }
        }
      }
    }
  }

  // Change cell.universe values from IDs to indices.
  for (Cell* c : cells) {
    auto search = universe_map.find(c->universe_);
    if (search != universe_map.end()) {
      c->universe_ = search->second;
    } else {
      std::stringstream err_msg;
      err_msg << "Could not find universe " << c->universe_
              << " specified on cell " << c->id_;
      fatal_error(err_msg);
    }
  }

  // Change all lattice universe values from IDs to indices.
  for (Lattice* l : lattices) {
    l->adjust_indices();
  }
}

//==============================================================================

void
assign_temperatures()
{
  for (Cell* c : cells) {
    // Ignore non-material cells and cells with defined temperature.
    if (c->material_.size() == 0) continue;
    if (c->sqrtkT_.size() > 0) continue;

    c->sqrtkT_.reserve(c->material_.size());
    for (auto i_mat : c->material_) {
      if (i_mat == MATERIAL_VOID) {
        // Set void region to 0K.
        c->sqrtkT_.push_back(0);

      } else {
        if (materials[i_mat]->temperature_ >= 0) {
          // This material has a default temperature; use that value.
          auto T = materials[i_mat]->temperature_;
          c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN * T));
        } else {
          // Use the global default temperature.
          c->sqrtkT_.push_back(std::sqrt(K_BOLTZMANN *
            settings::temperature_default));
        }
      }
    }
  }
}

//==============================================================================

int32_t
find_root_universe()
{
  // Find all the universes listed as a cell fill.
  std::unordered_set<int32_t> fill_univ_ids;
  for (Cell* c : cells) {
    fill_univ_ids.insert(c->fill_);
  }

  // Find all the universes contained in a lattice.
  for (Lattice* lat : lattices) {
    for (auto it = lat->begin(); it != lat->end(); ++it) {
      fill_univ_ids.insert(*it);
    }
    if (lat->outer_ != NO_OUTER_UNIVERSE) {
      fill_univ_ids.insert(lat->outer_);
    }
  }

  // Figure out which universe is not in the set.  This is the root universe.
  bool root_found {false};
  int32_t root_univ;
  for (int32_t i = 0; i < universes.size(); i++) {
    auto search = fill_univ_ids.find(universes[i]->id_);
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
neighbor_lists()
{
  write_message("Building neighboring cells lists for each surface...", 6);

  for (int i = 0; i < cells.size(); i++) {
    for (auto token : cells[i]->region_) {
      // Skip operator tokens.
      if (std::abs(token) >= OP_UNION) continue;

      // This token is a surface index.  Add the cell to the surface's list.
      if (token > 0) {
        surfaces[std::abs(token)-1]->neighbor_pos_.push_back(i);
      } else {
        surfaces[std::abs(token)-1]->neighbor_neg_.push_back(i);
      }
    }
  }

  for (Surface* surf : surfaces) {
    surf->neighbor_pos_.shrink_to_fit();
    surf->neighbor_neg_.shrink_to_fit();
  }
}

//==============================================================================

void
prepare_distribcell()
{
  // Find all cells listed in a DistribcellFilter.
  std::unordered_set<int32_t> distribcells;
  for (auto& filt : tally_filters) {
    auto* distrib_filt = dynamic_cast<DistribcellFilter*>(filt.get());
    if (distrib_filt) {
      distribcells.insert(distrib_filt->cell_);
    }
  }

  // Find all cells with distributed materials or temperatures.  Make sure that
  // the number of materials/temperatures matches the number of cell instances.
  for (int i = 0; i < cells.size(); i++) {
    Cell& c {*cells[i]};

    if (c.material_.size() > 1) {
      if (c.material_.size() != c.n_instances_) {
        std::stringstream err_msg;
        err_msg <<  "Cell " << c.id_ <<  " was specified with "
              << c.material_.size() << " materials but has " << c.n_instances_
              << " distributed instances. The number of materials must equal "
              "one or the number of instances.";
        fatal_error(err_msg);
      }
      distribcells.insert(i);
    }

    if (c.sqrtkT_.size() > 1) {
      if (c.sqrtkT_.size() != c.n_instances_) {
        std::stringstream err_msg;
        err_msg <<  "Cell " << c.id_ <<  " was specified with "
          << c.sqrtkT_.size() << " temperatures but has " << c.n_instances_
          << " distributed instances. The number of temperatures must equal "
          "one or the number of instances.";
        fatal_error(err_msg);
      }
      distribcells.insert(i);
    }
  }

  // Search through universes for distributed cells and assign each one a
  // unique distribcell array index.
  int distribcell_index = 0;
  std::vector<int32_t> target_univ_ids;
  for (Universe* u : universes) {
    for (auto cell_indx : u->cells_) {
      if (distribcells.find(cell_indx) != distribcells.end()) {
        cells[cell_indx]->distribcell_index_ = distribcell_index;
        target_univ_ids.push_back(u->id_);
        ++distribcell_index;
      }
    }
  }

  // Allocate the cell and lattice offset tables.
  int n_maps = target_univ_ids.size();
  for (Cell* c : cells) {
    if (c->type_ != FILL_MATERIAL) {
      c->offset_.resize(n_maps, C_NONE);
    }
  }
  for (Lattice* lat : lattices) {
    lat->allocate_offset_table(n_maps);
  }

  // Fill the cell and lattice offset tables.
  for (int map = 0; map < target_univ_ids.size(); map++) {
    auto target_univ_id = target_univ_ids[map];
    for (Universe* univ : universes) {
      int32_t offset {0};  // TODO: is this a bug?  It matches F90 implementation.
      for (int32_t cell_indx : univ->cells_) {
        Cell& c = *cells[cell_indx];

        if (c.type_ == FILL_UNIVERSE) {
          c.offset_[map] = offset;
          int32_t search_univ = c.fill_;
          offset += count_universe_instances(search_univ, target_univ_id);

        } else if (c.type_ == FILL_LATTICE) {
          Lattice& lat = *lattices[c.fill_];
          offset = lat.fill_offset_table(offset, target_univ_id, map);
        }
      }
    }
  }
}

//==============================================================================

void
count_cell_instances(int32_t univ_indx)
{
  for (int32_t cell_indx : universes[univ_indx]->cells_) {
    Cell& c = *cells[cell_indx];
    ++c.n_instances_;

    if (c.type_ == FILL_UNIVERSE) {
      // This cell contains another universe.  Recurse into that universe.
      count_cell_instances(c.fill_);

    } else if (c.type_ == FILL_LATTICE) {
      // This cell contains a lattice.  Recurse into the lattice universes.
      Lattice& lat = *lattices[c.fill_];
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
  if (universes[search_univ]->id_ == target_univ_id) {
    return 1;
  }

  int count {0};
  for (int32_t cell_indx : universes[search_univ]->cells_) {
    Cell& c = *cells[cell_indx];

    if (c.type_ == FILL_UNIVERSE) {
      int32_t next_univ = c.fill_;
      count += count_universe_instances(next_univ, target_univ_id);

    } else if (c.type_ == FILL_LATTICE) {
      Lattice& lat = *lattices[c.fill_];
      for (auto it = lat.begin(); it != lat.end(); ++it) {
        int32_t next_univ = *it;
        count += count_universe_instances(next_univ, target_univ_id);
      }
    }
  }

  return count;
}

//==============================================================================

std::string
distribcell_path_inner(int32_t target_cell, int32_t map, int32_t target_offset,
                       const Universe& search_univ, int32_t offset)
{
  std::stringstream path;

  path << "u" << search_univ.id_ << "->";

  // Check to see if this universe directly contains the target cell.  If so,
  // write to the path and return.
  for (int32_t cell_indx : search_univ.cells_) {
    if ((cell_indx == target_cell) && (offset == target_offset)) {
      Cell& c = *cells[cell_indx];
      path << "c" << c.id_;
      return path.str();
    }
  }

  // The target must be further down the geometry tree and contained in a fill
  // cell or lattice cell in this universe.  Find which cell contains the
  // target.
  std::vector<std::int32_t>::const_reverse_iterator cell_it
       {search_univ.cells_.crbegin()};
  for (; cell_it != search_univ.cells_.crend(); ++cell_it) {
    Cell& c = *cells[*cell_it];

    // Material cells don't contain other cells so ignore them.
    if (c.type_ != FILL_MATERIAL) {
      int32_t temp_offset;
      if (c.type_ == FILL_UNIVERSE) {
        temp_offset = offset + c.offset_[map];
      } else {
        Lattice& lat = *lattices[c.fill_];
        int32_t indx = lat.universes_.size()*map + lat.begin().indx_;
        temp_offset = offset + lat.offsets_[indx];
      }

      // The desired cell is the first cell that gives an offset smaller or
      // equal to the target offset.
      if (temp_offset <= target_offset) break;
    }
  }

  // Add the cell to the path string.
  Cell& c = *cells[*cell_it];
  path << "c" << c.id_ << "->";

  if (c.type_ == FILL_UNIVERSE) {
    // Recurse into the fill cell.
    offset += c.offset_[map];
    path << distribcell_path_inner(target_cell, map, target_offset,
                                   *universes[c.fill_], offset);
    return path.str();
  } else {
    // Recurse into the lattice cell.
    Lattice& lat = *lattices[c.fill_];
    path << "l" << lat.id_;
    for (ReverseLatticeIter it = lat.rbegin(); it != lat.rend(); ++it) {
      int32_t indx = lat.universes_.size()*map + it.indx_;
      int32_t temp_offset = offset + lat.offsets_[indx];
      if (temp_offset <= target_offset) {
        offset = temp_offset;
        path << "(" << lat.index_to_string(it.indx_) << ")->";
        path << distribcell_path_inner(target_cell, map, target_offset,
                                       *universes[*it], offset);
        return path.str();
      }
    }
    throw std::runtime_error{"Error determining distribcell path."};
  }
}

std::string
distribcell_path(int32_t target_cell, int32_t map, int32_t target_offset)
{
  auto& root_univ = *universes[openmc_root_universe];
  return distribcell_path_inner(target_cell, map, target_offset, root_univ, 0);
}

//==============================================================================

int
maximum_levels(int32_t univ)
{
  int levels_below {0};

  for (int32_t cell_indx : universes[univ]->cells_) {
    Cell& c = *cells[cell_indx];
    if (c.type_ == FILL_UNIVERSE) {
      int32_t next_univ = c.fill_;
      levels_below = std::max(levels_below, maximum_levels(next_univ));
    } else if (c.type_ == FILL_LATTICE) {
      Lattice& lat = *lattices[c.fill_];
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
  for (Cell* c : cells) {delete c;}
  cells.clear();
  cell_map.clear();
  n_cells = 0;

  for (Universe* u : universes) {delete u;}
  universes.clear();
  universe_map.clear();

  for (Lattice* lat : lattices) {delete lat;}
  lattices.clear();
  lattice_map.clear();

  overlap_check_count.clear();
}

} // namespace openmc
