#include "geometry.h"

#include <sstream>

#include "cell.h"
#include "constants.h"
#include "error.h"
#include "lattice.h"
#include "settings.h"


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
      //=======================================================================
      //! Found a material cell which means this is the lowest coord level.

      // Find the distribcell instance number.
      if (c.material.size() > 1 || c.sqrtkT.size() > 1) {
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

    } else if (c.type == FILL_UNIVERSE) {
      //========================================================================
      //! Found a lower universe, update this coord level then search the next.

      // Set the lower coordinate level universe.
      p->coord[p->n_coord].universe = c.fill + 1;

      // Set the position and direction.
      for (int i = 0; i < 3; i++) {
        p->coord[p->n_coord].xyz[i] = p->coord[p->n_coord-1].xyz[i];
        p->coord[p->n_coord].uvw[i] = p->coord[p->n_coord-1].uvw[i];
      }

      // Apply translation.
      p->coord[p->n_coord].xyz[0] -= c.translation.x;
      p->coord[p->n_coord].xyz[1] -= c.translation.y;
      p->coord[p->n_coord].xyz[2] -= c.translation.z;

      // Apply rotation.
      if (!c.rotation.empty()) {
        auto x = p->coord[p->n_coord].xyz[0];
        auto y = p->coord[p->n_coord].xyz[1];
        auto z = p->coord[p->n_coord].xyz[2];
        p->coord[p->n_coord].xyz[0] = x*c.rotation[3] + y*c.rotation[4]
                                      + z*c.rotation[5];
        p->coord[p->n_coord].xyz[1] = x*c.rotation[6] + y*c.rotation[7]
                                      + z*c.rotation[8];
        p->coord[p->n_coord].xyz[2] = x*c.rotation[9] + y*c.rotation[10]
                                      + z*c.rotation[11];
        auto u = p->coord[p->n_coord].uvw[0];
        auto v = p->coord[p->n_coord].uvw[1];
        auto w = p->coord[p->n_coord].uvw[2];
        p->coord[p->n_coord].uvw[0] = u*c.rotation[3] + v*c.rotation[4]
                                      + w*c.rotation[5];
        p->coord[p->n_coord].uvw[1] = u*c.rotation[6] + v*c.rotation[7]
                                      + w*c.rotation[8];
        p->coord[p->n_coord].uvw[2] = u*c.rotation[9] + v*c.rotation[10]
                                      + w*c.rotation[11];
        p->coord[p->n_coord].rotated = true;
      }

      // Update the coordinate level and recurse.
      ++p->n_coord;
      find_cell(p, 0, nullptr);

    } else if (c.type == FILL_LATTICE) {
      //========================================================================
      //! Found a lower lattice, update this coord level then search the next.

      Lattice& lat {*lattices_c[c.fill]};

      // Determine lattice indices.
      Position r {p->coord[p->n_coord-1].xyz};
      Direction u {p->coord[p->n_coord-1].uvw};
      r += TINY_BIT * u;
      auto i_xyz = lat.get_indices(r);

      // Store lower level coordinates.
      r = lat.get_local_position(p->coord[p->n_coord-1].xyz, i_xyz);
      p->coord[p->n_coord].xyz[0] = r.x;
      p->coord[p->n_coord].xyz[1] = r.y;
      p->coord[p->n_coord].xyz[2] = r.z;
      p->coord[p->n_coord].uvw[0] = u.x;
      p->coord[p->n_coord].uvw[1] = u.y;
      p->coord[p->n_coord].uvw[2] = u.z;

      // Set lattice indices.
      p->coord[p->n_coord].lattice = c.fill + 1;
      p->coord[p->n_coord].lattice_x = i_xyz[0];
      p->coord[p->n_coord].lattice_y = i_xyz[1];
      p->coord[p->n_coord].lattice_z = i_xyz[2];

      // Set the lower coordinate level universe.
      if (lat.are_valid_indices(i_xyz)) {
        p->coord[p->n_coord].universe = lat[i_xyz] + 1;
      } else {
        if (lat.outer != NO_OUTER_UNIVERSE) {
          p->coord[p->n_coord].universe = lat.outer + 1;
        } else {
          std::stringstream err_msg;
          err_msg << "Particle " << p->id << " is outside lattice "
                  << lat.id << " but the lattice has no defined outer "
                  "universe.";
          warning(err_msg);
          found = false;
        }
      }

      // Update the coordinate level and recurse.
      ++p->n_coord;
      find_cell(p, 0, nullptr);
    }
  }

  return found;
}

} // namespace openmc
