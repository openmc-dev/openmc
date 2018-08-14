#include "geometry.h"

#include <sstream>

#include "cell.h"
#include "error.h"
#include "particle.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

std::vector<int64_t> overlap_check_count;

extern "C" bool
check_cell_overlap(Particle* p) {
  int n_coord = p->n_coord;

  // loop through each coordinate level
  for (int j = 0; j < n_coord; j++) {
    //p->n_coord = j + 1;
    Universe& univ {*global_universes[p->coord[j].universe - 1]};
    int n = univ.cells.size();

    // loop through each cell on this level
    for (int i = 0; i < n; i++) {
      int index_cell = univ.cells[i];
      Cell& c {*global_cells[index_cell]};

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

} // namespace openmc
