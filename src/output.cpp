#include "openmc/output.h"

#include <algorithm>  // for std::transform
#include <cstring>  // for strlen
#include <ctime> // for time, localtime
#include <iomanip>  // for setw, setprecision, put_time
#include <ios> // for fixed, scientific, left
#include <iostream>
#include <sstream>
#include <utility> // for pair

#include <omp.h>
#include "xtensor/xview.hpp"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/plot.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"
#include "openmc/timer.h"
#include "openmc/tallies/tally.h"

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
    std::cout << '\n' << out.str() << "\n\n";
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
  // Save state of cout
  auto f {std::cout.flags()};

  // Determine overall generation and number of active generations
  int i = overall_generation() - 1;
  int n = simulation::current_batch > settings::n_inactive ?
    settings::gen_per_batch*n_realizations + simulation::current_gen : 0;

  // Set format for values
  std::cout << std::fixed << std::setprecision(5);

  // write out information batch and option independent output
  std::cout << "  "  << std::setw(9) <<  std::to_string(simulation::current_batch)
    + "/" + std::to_string(simulation::current_gen) << "   " << std::setw(8)
    << simulation::k_generation[i];

  // write out entropy info
  if (settings::entropy_on) {
    std::cout << "   " << std::setw(8) << simulation::entropy[i];
  }

  if (n > 1) {
    std::cout << "   " << std::setw(8) << simulation::keff << " +/-"
      << std::setw(8) << simulation::keff_std;
  }
  std::cout << '\n';

  // Restore state of cout
  std::cout.flags(f);
}

//==============================================================================

void print_batch_keff()
{
  // Save state of cout
  auto f {std::cout.flags()};

  // Determine overall generation and number of active generations
  int i = simulation::current_batch*settings::gen_per_batch - 1;
  int n = n_realizations*settings::gen_per_batch;

  // Set format for values
  std::cout << std::fixed << std::setprecision(5);

  // write out information batch and option independent output
  std::cout << "  "  << std::setw(9) <<  std::to_string(simulation::current_batch)
    + "/" + std::to_string(settings::gen_per_batch) << "   " << std::setw(8)
    << simulation::k_generation[i];

  // write out entropy info
  if (settings::entropy_on) {
    std::cout << "   " << std::setw(8) << simulation::entropy[i];
  }

  if (n > 1) {
    std::cout << "   " << std::setw(8) << simulation::keff << " +/-"
      << std::setw(8) << simulation::keff_std;
  }
  std::cout << '\n';

  // Restore state of cout
  std::cout.flags(f);
}

//==============================================================================

void show_time(const char* label, double secs, int indent_level=0)
{
  std::cout << std::string(2*indent_level, ' ');
  int width = 33 - indent_level*2;
  std::cout << " " << std::setw(width) << std::left << label << " = "
    << std::setw(10) << std::right << secs << " seconds\n";
}

void show_rate(const char* label, double particles_per_sec)
{
  std::cout << " " << std::setw(33) << std::left << label << " = " <<
    particles_per_sec << " particles/second\n";
}

void print_runtime()
{
  using namespace simulation;

  // display header block
  header("Timing Statistics", 6);

  // Save state of cout
  auto f {std::cout.flags()};

  // display time elapsed for various sections
  std::cout << std::scientific << std::setprecision(4);
  show_time("Total time for initialization", time_initialize.elapsed());
  show_time("Reading cross sections", time_read_xs.elapsed(), 1);
  show_time("Total time in simulation", time_inactive.elapsed() +
    time_active.elapsed());
  show_time("Time in transport only", time_transport.elapsed(), 1);
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    show_time("Time in inactive batches", time_inactive.elapsed(), 1);
  }
  show_time("Time in active batches", time_active.elapsed(), 1);
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    show_time("Time synchronizing fission bank", time_bank.elapsed(), 1);
    show_time("Sampling source sites", time_bank_sample.elapsed(), 2);
    show_time("SEND/RECV source sites", time_bank_sendrecv.elapsed(), 2);
  }
  show_time("Time accumulating tallies", time_tallies.elapsed(), 1);
  show_time("Total time for finalization", time_finalize.elapsed());
  show_time("Total time elapsed", time_total.elapsed());

  // Restore state of cout
  std::cout.flags(f);

  // Calculate particle rate in active/inactive batches
  int n_active = simulation::current_batch - settings::n_inactive;
  double speed_inactive;
  double speed_active;
  if (settings::restart_run) {
    if (simulation::restart_batch < settings::n_inactive) {
      speed_inactive = (settings::n_particles * (settings::n_inactive
        - simulation::restart_batch) * settings::gen_per_batch)
        / time_inactive.elapsed();
      speed_active = (settings::n_particles * n_active
        * settings::gen_per_batch) / time_active.elapsed();
    } else {
      speed_inactive = 0.0;
      speed_active = (settings::n_particles * (settings::n_batches
        - simulation::restart_batch) * settings::gen_per_batch)
        / time_active.elapsed();
    }
  } else {
    if (settings::n_inactive > 0) {
      speed_inactive = (settings::n_particles * settings::n_inactive
        * settings::gen_per_batch) / time_inactive.elapsed();
    }
    speed_active = (settings::n_particles * n_active * settings::gen_per_batch)
      / time_active.elapsed();
  }

  // display calculation rate
  std::cout << std::setprecision(6) << std::showpoint;
  if (!(settings::restart_run && (simulation::restart_batch >= settings::n_inactive))
      && settings::n_inactive > 0) {
    show_rate("Calculation Rate (inactive)", speed_inactive);
  }
  show_rate("Calculation Rate (active)", speed_active);

  // Restore state of cout
  std::cout.flags(f);
}

//==============================================================================

std::pair<double, double> mean_stdev(const double* x, int n)
{
  double mean = x[RESULT_SUM] / n;
  double stdev = n > 1 ? std::sqrt((x[RESULT_SUM_SQ]/n
    - mean*mean)/(n - 1)) : 0.0;
  return {mean, stdev};
}

//==============================================================================

void print_results()
{
  // Save state of cout
  auto f {std::cout.flags()};

  // display header block for results
  header("Results", 4);

  // Calculate t-value for confidence intervals
  int n = n_realizations;
  double alpha, t_n1, t_n3;
  if (settings::confidence_intervals) {
    alpha = 1.0 - CONFIDENCE_LEVEL;
    t_n1 = t_percentile(1.0 - alpha/2.0, n - 1);
    t_n3 = t_percentile(1.0 - alpha/2.0, n - 3);
  } else {
    t_n1 = 1.0;
    t_n3 = 1.0;
  }

  // Set formatting for floats
  std::cout << std::fixed << std::setprecision(5);

  // write global tallies
  auto gt = global_tallies();
  double mean, stdev;
  if (n > 1) {
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      std::tie(mean, stdev) = mean_stdev(&gt(K_COLLISION, 0), n);
      std::cout << " k-effective (Collision)     = "
        << mean << " +/- " << t_n1 * stdev << '\n';
      std::tie(mean, stdev) = mean_stdev(&gt(K_TRACKLENGTH, 0), n);
      std::cout << " k-effective (Track-length)  = "
        << mean << " +/- " << t_n1 * stdev << '\n';
      std::tie(mean, stdev) = mean_stdev(&gt(K_ABSORPTION, 0), n);
      std::cout << " k-effective (Absorption)    = "
        << mean << " +/- " << t_n1 * stdev << '\n';
      if (n > 3) {
        double k_combined[2];
        openmc_get_keff(k_combined);
        std::cout << " Combined k-effective        = "
          << k_combined[0] << " +/- " << t_n3 * k_combined[1] << '\n';
      }
    }
    std::tie(mean, stdev) = mean_stdev(&gt(LEAKAGE, 0), n);
    std::cout << " Leakage Fraction            = "
      << mean << " +/- " << t_n1 * stdev << '\n';
  } else {
    if (mpi::master) warning("Could not compute uncertainties -- only one "
      "active batch simulated!");

    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      std::cout << " k-effective (Collision)    = "
        << gt(K_COLLISION, RESULT_SUM) / n << '\n';
      std::cout << " k-effective (Track-length) = "
        << gt(K_TRACKLENGTH, RESULT_SUM) / n << '\n';
      std::cout << " k-effective (Absorption)   = "
        << gt(K_ABSORPTION, RESULT_SUM) / n << '\n';
    }
    std::cout << " Leakage Fraction           = "
      << gt(LEAKAGE, RESULT_SUM) / n << '\n';
  }
  std::cout << '\n';

  // Restore state of cout
  std::cout.flags(f);
}

} // namespace openmc
