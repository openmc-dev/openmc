#include "openmc/output.h"

#include <algorithm>  // for std::transform
#include <cstring>  // for strlen
#include <iomanip>  // for setw
#include <iostream>
#include <sstream>
#include <ctime>

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/capi.h"
#include "openmc/settings.h"
#include "openmc/plot.h"

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
  if (settings::verbosity >= level) {
    std::cout << out.str() << "\n\n";
  }
}

std::string time_stamp()
{
  int base_year = 1990;
  std::stringstream ts;
  std::time_t t = std::time(0);   // get time now
  std::tm* now = std::localtime(&t);
  ts << now->tm_year + base_year << "-" << now->tm_mon
     << "-" << now->tm_mday << " " << now->tm_hour
     << ":" << now->tm_min << ":" << now->tm_sec;
  return ts.str();
}

//==============================================================================

//===============================================================================
// PRINT_PLOT displays selected options for plotting
//===============================================================================

void print_plot() {

  header("PLOTTING SUMMARY", 5);

  for (auto pl : model::plots) {
    // Plot id
    std::cout << "Plot ID: " << pl.id_ << "\n";
    // Plot filename
    std::cout << "Plot file: " << pl.path_plot_ << "\n";
    // Plot level
    std::cout << "Universe depth: " << pl.level_ << "\n";

    // Plot type
    if (PlotType::slice == pl.type_) {
      std::cout << "Plot Type: Slice" << "\n";
    } else if (PlotType::voxel == pl.type_) {
      std::cout << "Plot Type: Voxel" << "\n";
    }

    // Plot parameters
    std::cout << "Origin: " << pl.origin_[0] << " "
              << pl.origin_[1] << " "
              << pl.origin_[2] << "\n";

    if (PlotType::slice == pl.type_) {
      std::cout << std::setprecision(4)
                << "Width: "
                << pl.width_[0] << " "
                << pl.width_[1] << "\n";
    } else if (PlotType::voxel == pl.type_) {
      std::cout << std::setprecision(4)
                << "Width: "
                << pl.width_[0] << " "
                << pl.width_[1] << " "
                << pl.width_[2] << "\n";
    }

    if (PlotColorBy::cells == pl.color_by_) {
      std::cout << "Coloring: Cells" << "\n";
    } else if (PlotColorBy::mats == pl.color_by_) {
      std::cout << "Coloring: Materials" << "\n";
    }

    if (PlotType::slice == pl.type_) {
      switch(pl.basis_) {
      case PlotBasis::xy:
        std::cout <<  "Basis: XY" << "\n";
        break;
      case PlotBasis::xz:
        std::cout <<  "Basis: XZ" << "\n";
        break;
      case PlotBasis::yz:
        std::cout <<  "Basis: YZ" << "\n";
        break;
      }
      std::cout << "Pixels: " << pl.pixels_[0] << " "
                << pl.pixels_[1] << " " << "\n";
    } else if (PlotType::voxel == pl.type_) {
      std::cout << "Voxels: " << pl.pixels_[0] << " "
                << pl.pixels_[1] << " "
                << pl.pixels_[2] << "\n";
    }

    std::cout << "\n";

  }
}

void
print_overlap_check() {
#ifdef OPENMC_MPI
  std::vector<int64_t> temp(model::overlap_check_count);
  int err = MPI_Reduce(temp.data(), model::overlap_check_count.data(),
                       model::overlap_check_count.size(), MPI_INT64_T, MPI_SUM, 0,
                       mpi::intracomm);
#endif

  if (mpi::master) {
    header("cell overlap check summary", 1);
    std::cout << " Cell ID      No. Overlap Checks\n";

    std::vector<int32_t> sparse_cell_ids;
    for (int i = 0; i < model::cells.size(); i++) {
      std::cout << " " << std::setw(8) << model::cells[i]->id_ << std::setw(17)
                << model::overlap_check_count[i] << "\n";
      if (model::overlap_check_count[i] < 10) {
        sparse_cell_ids.push_back(model::cells[i]->id_);
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
