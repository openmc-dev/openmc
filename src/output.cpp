#include "openmc/output.h"

#include <algorithm>  // for std::transform
#include <cstring>  // for strlen
#include <iomanip>  // for setw
#include <iostream>
#include <sstream>

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/capi.h"
#include "openmc/settings.h"


namespace openmc {

void
header(const char* msg, int level) {
  // Determine how many times to repeat the '=' character.
  int n_prefix = (63 - strlen(msg)) / 2;
  int n_suffix = n_prefix;
  if ((strlen(msg) % 2) == 0) ++n_suffix;

  // Convert to uppercase.
  std::string upper(msg);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

  // Add ===>  <=== markers.
  std::stringstream out;
  out << ' ';
  for (int i = 0; i < n_prefix; i++) out << '=';
  out << ">     " << upper << "     <";
  for (int i = 0; i < n_suffix; i++) out << '=';

  // Print header based on verbosity level.
  if (openmc_verbosity >= level) {
    std::cout << out.str() << std::endl << std::endl;
  }
}

//==============================================================================

extern "C" void
print_overlap_check() {
#ifdef OPENMC_MPI
  std::vector<int64_t> temp(overlap_check_count);
  int err = MPI_Reduce(temp.data(), overlap_check_count.data(),
                       overlap_check_count.size(), MPI_INT64_T, MPI_SUM, 0,
                       mpi::intracomm);
#endif

  if (openmc_master) {
    header("cell overlap check summary", 1);
    std::cout << " Cell ID      No. Overlap Checks\n";

    std::vector<int32_t> sparse_cell_ids;
    for (int i = 0; i < n_cells; i++) {
      std::cout << " " << std::setw(8) << global_cells[i]->id_ << std::setw(17)
                << overlap_check_count[i] << "\n";
      if (overlap_check_count[i] < 10) {
        sparse_cell_ids.push_back(global_cells[i]->id_);
      }
    }

    std::cout << "\n There were " << sparse_cell_ids.size()
              << " cells with less than 10 overlap checks\n";
    for (auto id : sparse_cell_ids) {
      std::cout << " " << id;
    }
    std::cout << "\n";
  }
}

} // namespace openmc
