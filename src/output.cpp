#include "openmc/output.h"

#include <algorithm>  // for transform, max
#include <cstring>  // for strlen
#include <ctime> // for time, localtime
#include <iomanip>  // for setw, setprecision, put_time
#include <ios> // for fixed, scientific, left
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility> // for pair

#include <fmt/core.h>
#include <fmt/ostream.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/plot.h"
#include "openmc/reaction.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/surface.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

namespace openmc {

//==============================================================================

void title()
{
  fmt::print(
    "                                %%%%%%%%%%%%%%%\n"
    "                           %%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                                    %%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                                     %%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                 ###############      %%%%%%%%%%%%%%%%%%%%%%%%\n"
    "                ##################     %%%%%%%%%%%%%%%%%%%%%%%\n"
    "                ###################     %%%%%%%%%%%%%%%%%%%%%%%\n"
    "                ####################     %%%%%%%%%%%%%%%%%%%%%%\n"
    "                #####################     %%%%%%%%%%%%%%%%%%%%%\n"
    "                ######################     %%%%%%%%%%%%%%%%%%%%\n"
    "                #######################     %%%%%%%%%%%%%%%%%%\n"
    "                 #######################     %%%%%%%%%%%%%%%%%\n"
    "                 ######################     %%%%%%%%%%%%%%%%%\n"
    "                  ####################     %%%%%%%%%%%%%%%%%\n"
    "                    #################     %%%%%%%%%%%%%%%%%\n"
    "                     ###############     %%%%%%%%%%%%%%%%\n"
    "                       ############     %%%%%%%%%%%%%%%\n"
    "                          ########     %%%%%%%%%%%%%%\n"
    "                                      %%%%%%%%%%%\n\n");

  // Write version information
  fmt::print(
    "                   | The OpenMC Monte Carlo Code\n"
    "         Copyright | 2011-2020 MIT and OpenMC contributors\n"
    "           License | http://openmc.readthedocs.io/en/latest/license.html\n"
    "           Version | {}.{}.{}{}\n", VERSION_MAJOR, VERSION_MINOR,
    VERSION_RELEASE, VERSION_DEV ? "-dev" : "");
#ifdef GIT_SHA1
  fmt::print("          Git SHA1 | {}\n", GIT_SHA1);
#endif

  // Write the date and time
  fmt::print("         Date/Time | {}\n", time_stamp());

#ifdef OPENMC_MPI
  // Write number of processors
  fmt::print("     MPI Processes | {}\n", mpi::n_procs);
#endif

#ifdef _OPENMP
  // Write number of OpenMP threads
  fmt::print("    OpenMP Threads | {}\n", omp_get_max_threads());
#endif
  std::cout << std::endl;
}

//==============================================================================

std::string
header(const char* msg) {
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

  return out.str();
}

std::string header(const std::string& msg) {return header(msg.c_str());}

void
header(const char* msg, int level) {
  auto out = header(msg);

  // Print header based on verbosity level.
  if (settings::verbosity >= level)
    std::cout << '\n' << out << "\n" << std::endl;
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
  switch (p->type_) {
    case Particle::Type::neutron:
      fmt::print("Neutron ");
      break;
    case Particle::Type::photon:
      fmt::print("Photon ");
      break;
    case Particle::Type::electron:
      fmt::print("Electron ");
      break;
    case Particle::Type::positron:
      fmt::print("Positron ");
      break;
    default:
      fmt::print("Unknown Particle ");
  }
  fmt::print("{}\n", p->id_);

  // Display particle geometry hierarchy.
  for (auto i = 0; i < p->n_coord_; i++) {
    fmt::print("  Level {}\n", i);

    if (p->coord_[i].cell != C_NONE) {
      const Cell& c {*model::cells[p->coord_[i].cell]};
      fmt::print("    Cell             = {}\n", c.id_);
    }

    if (p->coord_[i].universe != C_NONE) {
      const Universe& u {*model::universes[p->coord_[i].universe]};
      fmt::print("    Universe         = {}\n", u.id_);
    }

    if (p->coord_[i].lattice != C_NONE) {
      const Lattice& lat {*model::lattices[p->coord_[i].lattice]};
      fmt::print("    Lattice          = {}\n", lat.id_);
      fmt::print("    Lattice position = ({},{},{})\n", p->coord_[i].lattice_x,
        p->coord_[i].lattice_y, p->coord_[i].lattice_z);
    }

    fmt::print("    r = {}\n", p->coord_[i].r);
    fmt::print("    u = {}\n", p->coord_[i].u);
  }

  // Display miscellaneous info.
  if (p->surface_ != 0) {
    const Surface& surf {*model::surfaces[std::abs(p->surface_)-1]};
    fmt::print("  Surface = {}\n", (p->surface_ > 0) ? surf.id_ : -surf.id_);
  }
  fmt::print("  Weight = {}\n", p->wgt_);
  if (settings::run_CE) {
    fmt::print("  Energy = {}\n", p->E_);
  } else {
    fmt::print("  Energy Group = {}\n", p->g_);
  }
  fmt::print("  Delayed Group = {}\n\n", p->delayed_group_);
}

//==============================================================================

void print_plot()
{
  header("PLOTTING SUMMARY", 5);
  if (settings::verbosity < 5) return;

  for (auto pl : model::plots) {
    // Plot id
    fmt::print("Plot ID: {}\n", pl.id_);
    // Plot filename
    fmt::print("Plot file: {}\n", pl.path_plot_);
    // Plot level
    fmt::print("Universe depth: {}\n", pl.level_);

    // Plot type
    if (PlotType::slice == pl.type_) {
      fmt::print("Plot Type: Slice\n");
    } else if (PlotType::voxel == pl.type_) {
      fmt::print("Plot Type: Voxel\n");
    }

    // Plot parameters
    fmt::print("Origin: {} {} {}\n", pl.origin_[0], pl.origin_[1], pl.origin_[2]);

    if (PlotType::slice == pl.type_) {
      fmt::print("Width: {:4} {:4}\n", pl.width_[0], pl.width_[1]);
    } else if (PlotType::voxel == pl.type_) {
      fmt::print("Width: {:4} {:4} {:4}\n", pl.width_[0], pl.width_[1],
        pl.width_[2]);
    }

    if (PlotColorBy::cells == pl.color_by_) {
      fmt::print("Coloring: Cells\n");
    } else if (PlotColorBy::mats == pl.color_by_) {
      fmt::print("Coloring: Materials\n");
    }

    if (PlotType::slice == pl.type_) {
      switch(pl.basis_) {
      case PlotBasis::xy:
        fmt::print("Basis: XY\n");
        break;
      case PlotBasis::xz:
        fmt::print("Basis: XZ\n");
        break;
      case PlotBasis::yz:
        fmt::print("Basis: YZ\n");
        break;
      }
      fmt::print("Pixels: {} {}\n", pl.pixels_[0], pl.pixels_[1]);
    } else if (PlotType::voxel == pl.type_) {
      fmt::print("Voxels: {} {} {}\n", pl.pixels_[0], pl.pixels_[1], pl.pixels_[2]);
    }

    fmt::print("\n");
  }
}

//==============================================================================

void
print_overlap_check()
{
  #ifdef OPENMC_MPI
    std::vector<int64_t> temp(model::overlap_check_count);
    MPI_Reduce(temp.data(), model::overlap_check_count.data(),
               model::overlap_check_count.size(), MPI_INT64_T,
               MPI_SUM, 0, mpi::intracomm);
  #endif

  if (mpi::master) {
    header("cell overlap check summary", 1);
    fmt::print(" Cell ID      No. Overlap Checks\n");

    std::vector<int32_t> sparse_cell_ids;
    for (int i = 0; i < model::cells.size(); i++) {
      fmt::print(" {:8}{:17}\n", model::cells[i]->id_, model::overlap_check_count[i]);
      if (model::overlap_check_count[i] < 10) {
        sparse_cell_ids.push_back(model::cells[i]->id_);
      }
    }

    fmt::print("\n There were {} cells with less than 10 overlap checks\n",
      sparse_cell_ids.size());
    for (auto id : sparse_cell_ids) {
      fmt::print(" {}", id);
    }
    fmt::print("\n");
  }
}

//==============================================================================

void print_usage()
{
  if (mpi::master) {
    fmt::print(
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
      "  -e, --event            Run using event-based parallelism\n"
      "  -v, --version          Show version information\n"
      "  -h, --help             Show this message\n");
  }
}

//==============================================================================

void print_version()
{
  if (mpi::master) {
    fmt::print("OpenMC version {}.{}.{}\n", VERSION_MAJOR, VERSION_MINOR,
      VERSION_RELEASE);
#ifdef GIT_SHA1
    fmt::print("Git SHA1: {}\n", GIT_SHA1);
#endif
    fmt::print("Copyright (c) 2011-2019 Massachusetts Institute of "
      "Technology and OpenMC contributors\nMIT/X license at "
      "<http://openmc.readthedocs.io/en/latest/license.html>\n");
  }
}

//==============================================================================

void print_columns()
{
  if (settings::entropy_on) {
    fmt::print(
      "  Bat./Gen.      k       Entropy         Average k \n"
      "  =========   ========   ========   ====================\n");
  } else {
    fmt::print(
      "  Bat./Gen.      k            Average k\n"
      "  =========   ========   ====================\n");
  }
}

//==============================================================================

void print_generation()
{
  // Determine overall generation and number of active generations
  int i = overall_generation() - 1;
  int n = simulation::current_batch > settings::n_inactive ?
    settings::gen_per_batch*simulation::n_realizations + simulation::current_gen : 0;

  // write out information batch and option independent output
  auto batch_and_gen = std::to_string(simulation::current_batch) + "/" +
    std::to_string(simulation::current_gen);
  fmt::print("  {:>9}   {:8.5f}", batch_and_gen, simulation::k_generation[i]);

  // write out entropy info
  if (settings::entropy_on) {
    fmt::print("   {:8.5f}", simulation::entropy[i]);
  }

  if (n > 1) {
    fmt::print("   {:8.5f} +/-{:8.5f}", simulation::keff, simulation::keff_std);
  }
  std::cout << std::endl;
}

//==============================================================================

void print_batch_keff()
{
  // Determine overall generation and number of active generations
  int i = simulation::current_batch*settings::gen_per_batch - 1;
  int n = simulation::n_realizations*settings::gen_per_batch;

  // write out information batch and option independent output
  auto batch_and_gen = std::to_string(simulation::current_batch) + "/" +
    std::to_string(settings::gen_per_batch);
  fmt::print("  {:>9}   {:8.5f}", batch_and_gen, simulation::k_generation[i]);

  // write out entropy info
  if (settings::entropy_on) {
    fmt::print("   {:8.5f}", simulation::entropy[i]);
  }

  if (n > 1) {
    fmt::print("   {:8.5f} +/-{:8.5f}", simulation::keff, simulation::keff_std);
  }
  std::cout << std::endl;
}

//==============================================================================

void show_time(const char* label, double secs, int indent_level=0)
{
  int width = 33 - indent_level*2;
  fmt::print("{0:{1}} {2:<{3}} = {4:>10.4e} seconds\n",
    "", 2*indent_level, label, width, secs);
}

void show_rate(const char* label, double particles_per_sec)
{
  fmt::print(" {:<33} = {:.6} particles/second\n", label, particles_per_sec);
}

void print_runtime()
{
  using namespace simulation;

  // display header block
  header("Timing Statistics", 6);
  if (settings::verbosity < 6) return;

  // display time elapsed for various sections
  show_time("Total time for initialization", time_initialize.elapsed());
  show_time("Reading cross sections", time_read_xs.elapsed(), 1);
  show_time("Total time in simulation", time_inactive.elapsed() +
    time_active.elapsed());
  show_time("Time in transport only", time_transport.elapsed(), 1);
  if (settings::event_based) {
    show_time("Particle initialization", time_event_init.elapsed(), 2);
    show_time("XS lookups", time_event_calculate_xs.elapsed(), 2);
    show_time("Advancing", time_event_advance_particle.elapsed(), 2);
    show_time("Surface crossings", time_event_surface_crossing.elapsed(), 2);
    show_time("Collisions", time_event_collision.elapsed(), 2);
    show_time("Particle death", time_event_death.elapsed(), 2);
  }
  if (settings::run_mode == RunMode::EIGENVALUE) {
    show_time("Time in inactive batches", time_inactive.elapsed(), 1);
  }
  show_time("Time in active batches", time_active.elapsed(), 1);
  if (settings::run_mode == RunMode::EIGENVALUE) {
    show_time("Time synchronizing fission bank", time_bank.elapsed(), 1);
    show_time("Sampling source sites", time_bank_sample.elapsed(), 2);
    show_time("SEND/RECV source sites", time_bank_sendrecv.elapsed(), 2);
  }
  show_time("Time accumulating tallies", time_tallies.elapsed(), 1);
  show_time("Total time for finalization", time_finalize.elapsed());
  show_time("Total time elapsed", time_total.elapsed());

  // Calculate particle rate in active/inactive batches
  int n_active = simulation::current_batch - settings::n_inactive;
  double speed_inactive = 0.0;
  double speed_active;
  if (settings::restart_run) {
    if (simulation::restart_batch < settings::n_inactive) {
      speed_inactive = (settings::n_particles * (settings::n_inactive
        - simulation::restart_batch) * settings::gen_per_batch)
        / time_inactive.elapsed();
      speed_active = (settings::n_particles * n_active
        * settings::gen_per_batch) / time_active.elapsed();
    } else {
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
  if (!(settings::restart_run && (simulation::restart_batch >= settings::n_inactive))
      && settings::n_inactive > 0) {
    show_rate("Calculation Rate (inactive)", speed_inactive);
  }
  show_rate("Calculation Rate (active)", speed_active);
}

//==============================================================================

std::pair<double, double>
mean_stdev(const double* x, int n)
{
  double mean = x[static_cast<int>(TallyResult::SUM)] / n;
  double stdev = n > 1 ? std::sqrt(std::max(0.0, (
    x[static_cast<int>(TallyResult::SUM_SQ)]/n - mean*mean)/(n - 1))) : 0.0;
  return {mean, stdev};
}

//==============================================================================

void print_results()
{
  // display header block for results
  header("Results", 4);
  if (settings::verbosity < 4) return;

  // Calculate t-value for confidence intervals
  int n = simulation::n_realizations;
  double alpha, t_n1, t_n3;
  if (settings::confidence_intervals) {
    alpha = 1.0 - CONFIDENCE_LEVEL;
    t_n1 = t_percentile(1.0 - alpha/2.0, n - 1);
    t_n3 = t_percentile(1.0 - alpha/2.0, n - 3);
  } else {
    t_n1 = 1.0;
    t_n3 = 1.0;
  }

  // write global tallies
  const auto& gt = simulation::global_tallies;
  double mean, stdev;
  if (n > 1) {
    if (settings::run_mode == RunMode::EIGENVALUE) {
      std::tie(mean, stdev) = mean_stdev(&gt(GlobalTally::K_COLLISION, 0), n);
      fmt::print(" k-effective (Collision)     = {:.5f} +/- {:.5f}\n",
        mean, t_n1 * stdev);
      std::tie(mean, stdev) = mean_stdev(&gt(GlobalTally::K_TRACKLENGTH, 0), n);
      fmt::print(" k-effective (Track-length)  = {:.5f} +/- {:.5f}\n",
        mean, t_n1 * stdev);
      std::tie(mean, stdev) = mean_stdev(&gt(GlobalTally::K_ABSORPTION, 0), n);
      fmt::print(" k-effective (Absorption)    = {:.5f} +/- {:.5f}\n",
        mean, t_n1 * stdev);
      if (n > 3) {
        double k_combined[2];
        openmc_get_keff(k_combined);
        fmt::print(" Combined k-effective        = {:.5f} +/- {:.5f}\n",
          k_combined[0], k_combined[1]);
      }
    }
    std::tie(mean, stdev) = mean_stdev(&gt(GlobalTally::LEAKAGE, 0), n);
    fmt::print(" Leakage Fraction            = {:.5f} +/- {:.5f}\n",
      mean, t_n1 * stdev);
  } else {
    if (mpi::master) warning("Could not compute uncertainties -- only one "
      "active batch simulated!");

    if (settings::run_mode == RunMode::EIGENVALUE) {
      fmt::print(" k-effective (Collision)    = {:.5f}\n",
        gt(GlobalTally::K_COLLISION, TallyResult::SUM) / n);
      fmt::print(" k-effective (Track-length) = {:.5f}\n",
        gt(GlobalTally::K_TRACKLENGTH, TallyResult::SUM) / n);
      fmt::print(" k-effective (Absorption)   = {:.5f}\n",
        gt(GlobalTally::K_ABSORPTION, TallyResult::SUM) / n);
    }
    fmt::print(" Leakage Fraction           = {:.5f}\n",
      gt(GlobalTally::LEAKAGE, TallyResult::SUM) / n);
  }
  fmt::print("\n");
}

//==============================================================================

const std::unordered_map<int, const char*> score_names = {
  {SCORE_FLUX,               "Flux"},
  {SCORE_TOTAL,              "Total Reaction Rate"},
  {SCORE_SCATTER,            "Scattering Rate"},
  {SCORE_NU_SCATTER,         "Scattering Production Rate"},
  {SCORE_ABSORPTION,         "Absorption Rate"},
  {SCORE_FISSION,            "Fission Rate"},
  {SCORE_NU_FISSION,         "Nu-Fission Rate"},
  {SCORE_KAPPA_FISSION,      "Kappa-Fission Rate"},
  {SCORE_EVENTS,             "Events"},
  {SCORE_DECAY_RATE,         "Decay Rate"},
  {SCORE_DELAYED_NU_FISSION, "Delayed-Nu-Fission Rate"},
  {SCORE_PROMPT_NU_FISSION,  "Prompt-Nu-Fission Rate"},
  {SCORE_INVERSE_VELOCITY,   "Flux-Weighted Inverse Velocity"},
  {SCORE_FISS_Q_PROMPT,      "Prompt fission power"},
  {SCORE_FISS_Q_RECOV,       "Recoverable fission power"},
  {SCORE_CURRENT,            "Current"},
};

//! Create an ASCII output file showing all tally results.

void
write_tallies()
{
  if (model::tallies.empty()) return;

  // Open the tallies.out file.
  std::ofstream tallies_out;
  tallies_out.open("tallies.out", std::ios::out | std::ios::trunc);

  // Loop over each tally.
  for (auto i_tally = 0; i_tally < model::tallies.size(); ++i_tally) {
    const auto& tally {*model::tallies[i_tally]};

    // Write header block.
    std::string tally_header("TALLY " + std::to_string(tally.id_));
    if (!tally.name_.empty()) tally_header += ": " + tally.name_;
    fmt::print(tallies_out, "{}\n\n", header(tally_header));

    if (!tally.writable_) {
      fmt::print(tallies_out, " Internal\n\n");
      continue;
    }

    // Calculate t-value for confidence intervals
    double t_value = 1;
    if (settings::confidence_intervals) {
      auto alpha = 1 - CONFIDENCE_LEVEL;
      t_value = t_percentile(1 - alpha*0.5, tally.n_realizations_ - 1);
    }

    // Write derivative information.
    if (tally.deriv_ != C_NONE) {
      const auto& deriv {model::tally_derivs[tally.deriv_]};
      switch (deriv.variable) {
      case DerivativeVariable::DENSITY:
        fmt::print(tallies_out, " Density derivative Material {}\n",
          deriv.diff_material);
        break;
      case DerivativeVariable::NUCLIDE_DENSITY:
        fmt::print(tallies_out, " Nuclide density derivative Material {} Nuclide {}\n",
          deriv.diff_material, data::nuclides[deriv.diff_nuclide]->name_);
        break;
      case DerivativeVariable::TEMPERATURE:
        fmt::print(tallies_out, " Temperature derivative Material {}\n",
          deriv.diff_material);
        break;
      default:
        fatal_error(fmt::format("Differential tally dependent variable for "
          "tally {} not defined in output.cpp", tally.id_));
      }
    }

    // Initialize Filter Matches Object
    std::vector<FilterMatch> filter_matches;
    // Allocate space for tally filter matches
    filter_matches.resize(model::tally_filters.size());

    // Loop over all filter bin combinations.
    auto filter_iter = FilterBinIter(tally, false, &filter_matches);
    auto end = FilterBinIter(tally, true, &filter_matches);
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;

      // Print info about this combination of filter bins.  The stride check
      // prevents redundant output.
      int indent = 0;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        if (filter_index % tally.strides(i) == 0) {
          auto i_filt = tally.filters(i);
          const auto& filt {*model::tally_filters[i_filt]};
          auto& match {filter_matches[i_filt]};
          fmt::print(tallies_out, "{0:{1}}{2}\n", "", indent + 1,
            filt.text_label(match.i_bin_));
        }
        indent += 2;
      }

      // Loop over all nuclide and score combinations.
      int score_index = 0;
      for (auto i_nuclide : tally.nuclides_) {
        // Write label for this nuclide bin.
        if (i_nuclide == -1) {
          fmt::print(tallies_out, "{0:{1}}Total Material\n", "", indent + 1);
        } else {
          if (settings::run_CE) {
            fmt::print(tallies_out, "{0:{1}}{2}\n", "", indent + 1,
              data::nuclides[i_nuclide]->name_);
          } else {
            fmt::print(tallies_out, "{0:{1}}{2}\n", "", indent + 1,
              data::mg.nuclides_[i_nuclide].name);
          }
        }

        // Write the score, mean, and uncertainty.
        indent += 2;
        for (auto score : tally.scores_) {
          std::string score_name = score > 0 ? reaction_name(score)
            : score_names.at(score);
          double mean, stdev;
          std::tie(mean, stdev) = mean_stdev(
            &tally.results_(filter_index, score_index, 0), tally.n_realizations_);
          fmt::print(tallies_out, "{0:{1}}{2:<36} {3:.6} +/- {4:.6}\n",
            "", indent + 1, score_name, mean, t_value * stdev);
          score_index += 1;
        }
        indent -= 2;
      }
    }
  }
}

} // namespace openmc
