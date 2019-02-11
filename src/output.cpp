#include "openmc/output.h"

#include <algorithm>  // for std::transform
#include <cstring>  // for strlen
#include <ctime>
#include <iomanip>  // for setw, setprecision
#include <ios>  // for left
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
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
#include "openmc/surface.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

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
    std::cout << out << "\n\n";
}

//==============================================================================

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

std::pair<double, double>
mean_stdev(double sum, double sum_sq, int n)
{
  double mean, std_dev;
  mean = sum / n;
  if (n > 1) {
    std_dev = std::sqrt((sum_sq / n - mean*mean) / (n - 1));
  } else {
    std_dev = 0;
  }
  return {mean, std_dev};
}

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

extern "C" void
write_tallies()
{
  if (model::tallies.empty()) return;

  // Open the tallies.out file.
  std::ofstream tallies_out;
  tallies_out.open("tallies.out", std::ios::out | std::ios::trunc);
  tallies_out << std::setprecision(6);

  // Loop over each tally.
  for (auto i_tally = 0; i_tally < model::tallies.size(); ++i_tally) {
    const auto& tally {*model::tallies[i_tally]};
    auto results = tally_results(i_tally+1);
    // TODO: get this directly from the tally object when it's been translated
    int32_t n_realizations;
    auto err = openmc_tally_get_n_realizations(i_tally+1, &n_realizations);

    // Calculate t-value for confidence intervals
    double t_value = 1;
    if (settings::confidence_intervals) {
      auto alpha = 1 - CONFIDENCE_LEVEL;
      t_value = t_percentile_c(1 - alpha*0.5, n_realizations - 1);
    }

    // Write header block.
    std::string tally_header("TALLY " + std::to_string(tally.id_));
    if (!tally.name_.empty()) tally_header += ": " + tally.name_;
    tallies_out << "\n" << header(tally_header) << "\n\n";

    // Write derivative information.
    if (tally.deriv_ != C_NONE) {
      const auto& deriv {model::tally_derivs[tally.deriv_]};
      switch (deriv.variable) {
      case DIFF_DENSITY:
        tallies_out << " Density derivative  Material "
          << std::to_string(deriv.diff_material) << "\n";
        break;
      case DIFF_NUCLIDE_DENSITY:
        tallies_out << " Nuclide density derivative  Material "
          << std::to_string(deriv.diff_material) << "  Nuclide "
          // TODO: off-by-one
          << data::nuclides[deriv.diff_nuclide-1]->name_ << "\n";
        break;
      case DIFF_TEMPERATURE:
        tallies_out << " Temperature derivative  Material "
          << std::to_string(deriv.diff_material) << "\n";
        break;
      default:
        fatal_error("Differential tally dependent variable for tally "
          + std::to_string(tally.id_) + " not defined in output.cpp");
      }
    }

    // Loop over all filter bin combinations.
    auto filter_iter = FilterBinIter(tally, false);
    auto end = FilterBinIter(tally, true);
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;

      // Print info about this combination of filter bins.  The stride check
      // prevents redundant output.
      int indent = 0;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        if ((filter_index-1) % tally.strides(i) == 0) {
          auto i_filt = tally.filters(i);
          const auto& filt {*model::tally_filters[i_filt]};
          auto& match {simulation::filter_matches[i_filt]};
          tallies_out << std::string(indent+1, ' ')
            << filt.text_label(match.i_bin_) << "\n";
        }
        indent += 2;
      }

      // Loop over all nuclide and score combinations.
      int score_index = 0;
      for (auto i_nuclide : tally.nuclides_) {
        // Write label for this nuclide bin.
        if (i_nuclide == -1) {
          tallies_out << std::string(indent+1, ' ') << "Total Material\n";
        } else {
          if (settings::run_CE) {
            tallies_out << std::string(indent+1, ' ')
              << data::nuclides[i_nuclide]->name_ << "\n";
          } else {
            tallies_out << std::string(indent+1, ' ')
              << data::nuclides_MG[i_nuclide].name << "\n";
          }
        }

        // Write the score, mean, and uncertainty.
        indent += 2;
        for (auto score : tally.scores_) {
          std::string score_name = score > 0 ? reaction_name(score)
            : score_names.at(score);
          double mean, stdev;
          //TODO: off-by-one
          std::tie(mean, stdev) = mean_stdev(
            results(filter_index-1, score_index, RESULT_SUM),
            results(filter_index-1, score_index, RESULT_SUM_SQ),
            n_realizations);
          tallies_out << std::string(indent+1, ' ')  << std::left
            << std::setw(36) << score_name << " " << mean << " +/- "
            << t_value * stdev << "\n";
          score_index += 1;
        }
        indent -= 2;
      }
    }
  }
}

} // namespace openmc
