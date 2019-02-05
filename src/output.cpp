#include "openmc/output.h"

#include <algorithm>  // for std::transform
#include <cstring>  // for strlen
#include <iomanip>  // for setw
#include <iostream>
#include <sstream>
#include <ctime>

#include <omp.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/message_passing.h"
#include "openmc/plot.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"

namespace openmc {

//==============================================================================

void title()
{
  std::cout <<
    "                                %%%%%%%%%%%%%%%\n" <<
    "                           %%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                                    %%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                                     %%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                 ###############      %%%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                ##################     %%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                ###################     %%%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                ####################     %%%%%%%%%%%%%%%%%%%%%%\n" <<
    "                #####################     %%%%%%%%%%%%%%%%%%%%%\n" <<
    "                ######################     %%%%%%%%%%%%%%%%%%%%\n" <<
    "                #######################     %%%%%%%%%%%%%%%%%%\n" <<
    "                 #######################     %%%%%%%%%%%%%%%%%\n" <<
    "                 ######################     %%%%%%%%%%%%%%%%%\n" <<
    "                  ####################     %%%%%%%%%%%%%%%%%\n" <<
    "                    #################     %%%%%%%%%%%%%%%%%\n" <<
    "                     ###############     %%%%%%%%%%%%%%%%\n" <<
    "                       ############     %%%%%%%%%%%%%%%\n" <<
    "                          ########     %%%%%%%%%%%%%%\n" <<
    "                                      %%%%%%%%%%%\n";

  // Write version information
  std::cout <<
    "                   | The OpenMC Monte Carlo Code\n" <<
    "         Copyright | 2011-2019 MIT and OpenMC contributors\n" <<
    "           License | http://openmc.readthedocs.io/en/latest/license.html\n" <<
    "           Version | " << VERSION_MAJOR << '.' << VERSION_MINOR << '.'
    << VERSION_RELEASE << '\n';
#ifdef GIT_SHA1
  std::cout << "          Git SHA1 | " << GIT_SHA1 << '\n';
#endif

  // Write the date and time
  std::cout << "         Date/Time | " << time_stamp() << '\n';

#ifdef OPENMC_MPI
  // Write number of processors
  std::cout << "     MPI Processes | " << mpi::n_procs << '\n';
#endif

#ifdef _OPENMP
  // Write number of OpenMP threads
  std::cout << "    OpenMC Threads | " << omp_get_max_threads() << '\n';
#endif
  std::cout << '\n';
}

//==============================================================================

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

//==============================================================================

std::string time_stamp()
{
  std::stringstream ts;
  std::time_t t = std::time(nullptr);   // get time now
  ts << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S");
  return ts.str();
}

//==============================================================================

extern "C" void print_particle(Particle* p)
{
  // Display particle type and ID.
  switch (p->type) {
    case static_cast<int>(ParticleType::neutron):
      std::cout << "Neutron ";
      break;
    case static_cast<int>(ParticleType::photon):
      std::cout << "Photon ";
      break;
    case static_cast<int>(ParticleType::electron):
      std::cout << "Electron ";
      break;
    case static_cast<int>(ParticleType::positron):
      std::cout << "Positron ";
      break;
    default:
      std::cout << "Unknown Particle ";
  }
  std::cout << p->id << "\n";

  // Display particle geometry hierarchy.
  for (auto i = 0; i < p->n_coord; i++) {
    std::cout << "  Level " << i << "\n";

    if (p->coord[i].cell != C_NONE) {
      const Cell& c {*model::cells[p->coord[i].cell]};
      std::cout << "    Cell             = " << c.id_ << "\n";
    }

    if (p->coord[i].universe != C_NONE) {
      const Universe& u {*model::universes[p->coord[i].universe]};
      std::cout << "    Universe         = " << u.id_ << "\n";
    }

    if (p->coord[i].lattice != F90_NONE) {
      const Lattice& lat {*model::lattices[p->coord[i].lattice]};
      std::cout << "    Lattice          = " << lat.id_ << "\n";
      std::cout << "    Lattice position = (" << p->coord[i].lattice_x
                << "," << p->coord[i].lattice_y << ","
                << p->coord[i].lattice_z << ")\n";
    }

    std::cout << "    xyz = " << p->coord[i].xyz[0] << " "
              << p->coord[i].xyz[1] << " " << p->coord[i].xyz[2] << "\n";
    std::cout << "    uvw = " << p->coord[i].uvw[0] << " "
              << p->coord[i].uvw[1] << " " << p->coord[i].uvw[2] << "\n";
  }

  // Display miscellaneous info.
  if (p->surface != ERROR_INT) {
    const Surface& surf {*model::surfaces[std::abs(p->surface)-1]};
    std::cout << "  Surface = " << std::copysign(surf.id_, p->surface) << "\n";
  }
  std::cout << "  Weight = " << p->wgt << "\n";
  if (settings::run_CE) {
    std::cout << "  Energy = " << p->E << "\n";
  } else {
    std::cout << "  Energy Group = " << p->g << "\n";
  }
  std::cout << "  Delayed Group = " << p->delayed_group << "\n";

  std::cout << "\n";
}

//==============================================================================

void print_plot()
{
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

//==============================================================================

void
print_overlap_check()
{
  #ifdef OPENMC_MPI
    std::vector<int64_t> temp(model::overlap_check_count);
    int err = MPI_Reduce(temp.data(), model::overlap_check_count.data(),
                         model::overlap_check_count.size(), MPI_INT64_T,
                         MPI_SUM, 0, mpi::intracomm);
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

//==============================================================================

void print_usage()
{
  if (mpi::master) {
    std::cout <<
      "Usage: openmc [options] [directory]\n\n"
      "Options:\n"
      "  -c, --volume           Run in stochastic volume calculation mode\n"
      "  -g, --geometry-debug   Run with geometry debugging on\n"
      "  -n, --particles        Number of particles per generation\n"
      "  -p, --plot             Run in plotting mode\n"
      "  -r, --restart          Restart a previous run from a state point\n"
      "                         or a particle restart file\n"
      "  -s, --threads          Number of OpenMP threads\n"
      "  -t, --track            Write tracks for all particles\n"
      "  -v, --version          Show version information\n"
      "  -h, --help             Show this message\n";
  }
}

//==============================================================================

void print_version()
{
  if (mpi::master) {
    std::cout << "OpenMC version " << VERSION_MAJOR << '.' << VERSION_MINOR
      << '.' << VERSION_RELEASE << '\n';
#ifdef GIT_SHA1
    std::cout << "Git SHA1: " << GIT_SHA1 << '\n';
#endif
    std::cout << "Copyright (c) 2011-2019 Massachusetts Institute of "
      "Technology and OpenMC contributors\nMIT/X license at "
      "<http://openmc.readthedocs.io/en/latest/license.html>\n";
  }
}

//==============================================================================

void print_columns()
{
  if (settings::entropy_on) {
    std::cout <<
      "  Bat./Gen.      k       Entropy         Average k \n"
      "  =========   ========   ========   ====================\n";
  } else {
    std::cout <<
      "  Bat./Gen.      k            Average k\n"
      "  =========   ========   ====================\n";
  }
}

//==============================================================================

void print_generation()
{
  // Determine overall generation and number of active generations
  int i = overall_generation() - 1;
  int n = simulation::current_batch > settings::n_inactive ?
    settings::gen_per_batch*n_realizations + simulation::current_gen : 0;

  // write out information batch and option independent output
  std::cout << "  "  << std::setw(9) <<  std::to_string(simulation::current_batch)
    + "/" + std::to_string(simulation::current_gen) << "   " <<
    std::fixed << std::setw(8) << std::setprecision(5) << simulation::k_generation[i];

  // write out entropy info
  if (settings::entropy_on) std::cout << "   " <<
    std::fixed << std::setw(8) << std::setprecision(5) << simulation::entropy[i];

  if (n > 1) {
    std::cout << "   " <<
      std::fixed << std::setw(8) << std::setprecision(5) << simulation::keff << " +/-" <<
      std::fixed << std::setw(8) << std::setprecision(5) << simulation::keff_std;
  }
  std::cout << '\n';
}

//==============================================================================

void print_batch_keff()
{
  // Determine overall generation and number of active generations
  int i = simulation::current_batch*settings::gen_per_batch - 1;
  int n = n_realizations*settings::gen_per_batch;

  // write out information batch and option independent output
  std::cout << "  "  << std::setw(9) <<  std::to_string(simulation::current_batch)
    + "/" + std::to_string(settings::gen_per_batch) << "   " <<
    std::fixed << std::setw(8) << std::setprecision(5) << simulation::k_generation[i];

  // write out entropy info
  if (settings::entropy_on) std::cout << "   " <<
    std::fixed << std::setw(8) << std::setprecision(5) << simulation::entropy[i];

  if (n > 1) {
    std::cout << "   " <<
      std::fixed << std::setw(8) << std::setprecision(5) << simulation::keff << " +/-" <<
      std::fixed << std::setw(8) << std::setprecision(5) << simulation::keff_std;
  }
  std::cout << '\n';
}

} // namespace openmc
