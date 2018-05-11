#include <sstream>
#include <unordered_set>

#include "cell.h"
#include "constants.h"
#include "error.h"
#include "lattice.h"

#include <iostream> //TODO: remove this


namespace openmc {

extern "C" int32_t
find_root_universe()
{
  // Find all the universes listed as a cell fill.
  std::unordered_set<int32_t> fill_univ_ids;
  for (Cell *c : cells_c) {
    fill_univ_ids.insert(c->fill);
  }

  // Find all the universes contained in a lattice.
  for (Lattice *lat : lattices_c) {
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
  for (int32_t i = 0; i < universes_c.size(); i++) {
    auto search = fill_univ_ids.find(universes_c[i]->id);
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

extern "C" void
adjust_indices_c()
{
  // Adjust material/fill idices.
  for (Cell *c : cells_c) {
    if (c->material[0] == C_NONE) {
      int32_t id = c->fill;
      auto search_univ = universe_dict.find(id);
      auto search_lat = lattice_dict.find(id);
      if (search_univ != universe_dict.end()) {
        c->type = FILL_UNIVERSE;
        c->fill = search_univ->second + 1;  //TODO: off-by-one
      } else if (search_lat != lattice_dict.end()) {
        c->type = FILL_LATTICE;
        c->fill = search_lat->second + 1;  //TODO: off-by-one
      } else {
        std::stringstream err_msg;
        err_msg << "Specified fill " << id << " on cell " << c->id
                << " is neither a universe nor a lattice.";
        fatal_error(err_msg);
      }
    } else {
      //TODO: materials
      c->type = FILL_MATERIAL;
    }
  }

  // Change cell.universe values from IDs to indices.
  for (Cell *c : cells_c) {
    auto search = universe_dict.find(c->universe);
    if (search != universe_dict.end()) {
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
  for (Lattice *l : lattices_c) {
    l->adjust_indices();
  }
}

} // namespace openmc
