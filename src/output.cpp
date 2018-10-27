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
    std::cout << out.str() << std::endl << std::endl;
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

  for (auto pl : plots) {
    // Plot id
    std::cout << "Plot ID: " << pl->id << std::endl;
    // Plot filename
    std::cout << "Plot file: " << pl->path_plot << std::endl;
    // Plot level
    std::cout << "Universe depth: " << pl->level << std::endl;

    // Plot type
    if (PLOT_TYPE::SLICE == pl->type) {
      std::cout << "Plot Type: Slice" << std::endl;
    } else if (PLOT_TYPE::VOXEL == pl->type) {
      std::cout << "Plot Type: Voxel" << std::endl;      
    }

    // Plot parameters
    std::cout << "Origin: " << pl->origin[0] << " "
              << pl->origin[1] << " "
              << pl->origin[2] << std::endl;

    if (PLOT_TYPE::SLICE == pl->type) {
      std::cout << std::setprecision(4)
                << "Width: "
                << pl->width[0] << " "
                << pl->width[1] << std::endl;
    } else if (PLOT_TYPE::VOXEL == pl->type) {
      std::cout << std::setprecision(4)
                << "Width: "
                << pl->width[0] << " "
                << pl->width[1] << " "
                << pl->width[2] << std::endl;
    }

    if (PLOT_COLOR_BY::CELLS == pl->color_by) {
      std::cout << "Coloring: Cells" << std::endl;
    } else if (PLOT_COLOR_BY::MATS == pl->color_by) {
      std::cout << "Coloring: Materials" << std::endl;      
    }
    
    if (PLOT_TYPE::SLICE == pl->type) {
      switch(pl->basis) {
      case PLOT_BASIS::XY:
        std::cout <<  "Basis: XY" << std::endl;
        break;
      case PLOT_BASIS::XZ:
        std::cout <<  "Basis: XZ" << std::endl;
        break;
      case PLOT_BASIS::YZ:
        std::cout <<  "Basis: YZ" << std::endl;
        break;
      }
      std::cout << "Pixels: " << pl->pixels[0] << " "
                << pl->pixels[1] << " " << std::endl;
    } else if (PLOT_TYPE::VOXEL == pl->type) {
      std::cout << "Voxels: " << pl->pixels[0] << " "
                << pl->pixels[1] << " "
                << pl->pixels[2] << std::endl;
    }

    std::cout << std::endl;
    
  }
}
  
void
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
      std::cout << " " << std::setw(8) << cells[i]->id_ << std::setw(17)
                << overlap_check_count[i] << "\n";
      if (overlap_check_count[i] < 10) {
        sparse_cell_ids.push_back(cells[i]->id_);
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
