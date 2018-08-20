#include "geometry.h"

#include <sstream>

#include "cell.h"
#include "constants.h"
#include "error.h"
#include "lattice.h"
#include "settings.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

std::vector<int64_t> overlap_check_count;

//==============================================================================

extern "C" bool
check_cell_overlap(Particle* p) {
  int n_coord = p->n_coord;

  // loop through each coordinate level
  for (int j = 0; j < n_coord; j++) {
    Universe& univ = *global_universes[p->coord[j].universe - 1];
    int n = univ.cells.size();

    // loop through each cell on this level
    for (int i = 0; i < n; i++) {
      int index_cell = univ.cells[i];
      Cell& c = *global_cells[index_cell];

      if (c.contains(p->coord[j].xyz, p->coord[j].uvw, p->surface)) {
        //TODO: off-by-one indexing
        if (index_cell != p->coord[j].cell - 1) {
          std::stringstream err_msg;
          err_msg << "Overlapping cells detected: " << c.id << ", "
                  << global_cells[p->coord[j].cell-1]->id << " on universe "
                  << univ.id;
          fatal_error(err_msg);
        }
        ++overlap_check_count[index_cell];
      }
    }
  }

  return false;
}

//==============================================================================

extern "C" bool
find_cell(Particle* p, int n_search_cells, int* search_cells) {
  for (int i = p->n_coord; i < MAX_COORD; i++) {
    p->coord[i].reset();
  }

  // Determine universe (if not yet set, use root universe)
  int i_universe = p->coord[p->n_coord-1].universe;
  if (i_universe == C_NONE) {
    p->coord[p->n_coord-1].universe = openmc_root_universe;
    i_universe = openmc_root_universe;
  }
  //TODO: off-by-one indexing
  --i_universe;

  // If not given a set of search cells, search all cells in the uninverse.
  if (n_search_cells == 0) {
    search_cells = global_universes[i_universe]->cells.data();
    n_search_cells = global_universes[i_universe]->cells.size();
  }

  // Find which cell of this universe the particle is in.
  bool found = false;
  int32_t i_cell;
  for (int i = 0; i < n_search_cells; i++) {
    i_cell = search_cells[i];

    // Make sure the search cell is in the same universe.
    //TODO: off-by-one indexing
    if (global_cells[i_cell]->universe - 1 != i_universe) continue;

    Position r {p->coord[p->n_coord-1].xyz};
    Direction u {p->coord[p->n_coord-1].uvw};
    int32_t surf = p->surface;
    if (global_cells[i_cell]->contains(r, u, surf)) {
      //TODO: off-by-one indexing
      p->coord[p->n_coord-1].cell = i_cell + 1;

      if (openmc_verbosity >= 10 || openmc_trace) {
        std::stringstream msg;
        msg << "    Entering cell " << global_cells[i_cell]->id;
        write_message(msg, 1);
      }
      found = true;
      break;
    }
  }

  if (found) {
    Cell& c {*global_cells[i_cell]};
    if (c.type == FILL_MATERIAL) {
      // Find the distribcell instance number.
      if (c.material.size() > 1 || c.sqrtkT.size() > 1) {
        //=====================================================================
        //! Found a material cell which means this is the lowest coord level.

        //TODO: off-by-one indexing
        int distribcell_index = c.distribcell_index - 1;
        int offset = 0;
        for (int i = 0; i < p->n_coord; i++) {
          Cell& c_i {*global_cells[p->coord[i].cell-1]};
          if (c_i.type == FILL_UNIVERSE) {
            offset += c_i.offset[distribcell_index];
          } else if (c_i.type == FILL_LATTICE) {
            Lattice& lat {*lattices_c[p->coord[i+1].lattice-1]};
            int i_xyz[3] {p->coord[i+1].lattice_x,
                          p->coord[i+1].lattice_y,
                          p->coord[i+1].lattice_z};
            if (lat.are_valid_indices(i_xyz)) {
              offset += lat.offset(distribcell_index, i_xyz);
            }
          }
        }
        p->cell_instance = offset + 1;
      } else {
        p->cell_instance = 1;
      }

      // Set the material and temperature.
      p->last_material = p->material;
      int32_t mat;
      if (c.material.size() > 1) {
        mat = c.material[p->cell_instance-1];
      } else {
        mat = c.material[0];
      }
      if (mat == MATERIAL_VOID) {
        p->material = MATERIAL_VOID;
      } else {
        p->material = mat + 1;
      }
      p->last_sqrtkT = p->sqrtkT;
      if (c.sqrtkT.size() > 1) {
        p->sqrtkT = c.sqrtkT[p->cell_instance-1];
      } else {
        p->sqrtkT = c.sqrtkT[0];
      }
    }
  }

  return found;
}

} // namespace openmc
