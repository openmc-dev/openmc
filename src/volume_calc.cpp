#include "openmc/volume_calc.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/timer.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"

#include <algorithm> // for copy
#include <cmath> // for pow, sqrt
#include <unordered_set>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<VolumeCalculation> volume_calcs;
}

//==============================================================================
// VolumeCalculation implementation
//==============================================================================

VolumeCalculation::VolumeCalculation(pugi::xml_node node)
{
  // Read domain type (cell, material or universe)
  std::string domain_type = get_node_value(node, "domain_type");
  if (domain_type == "cell") {
    domain_type_ = TallyDomain::CELL;
  } else if (domain_type == "material") {
    domain_type_ = TallyDomain::MATERIAL;
  } else if (domain_type == "universe") {
    domain_type_ = TallyDomain::UNIVERSE;
  } else {
    fatal_error(std::string("Unrecognized domain type for stochastic "
      "volume calculation: " + domain_type));
  }

  // Read domain IDs, bounding corodinates and number of samples
  domain_ids_ = get_node_array<int>(node, "domain_ids");
  lower_left_ = get_node_array<double>(node, "lower_left");
  upper_right_ = get_node_array<double>(node, "upper_right");
  n_samples_ = std::stoull(get_node_value(node, "samples"));

  if (check_for_node(node, "threshold")) {
    pugi::xml_node threshold_node = node.child("threshold");

    threshold_ = std::stod(get_node_value(threshold_node, "threshold"));
    if (threshold_ <= 0.0) {
      fatal_error(fmt::format("Invalid error threshold {} provided for a "
        "volume calculation.", threshold_));
    }

    std::string tmp = get_node_value(threshold_node, "type");
    if (tmp == "variance") {
      trigger_type_ = TriggerMetric::variance;
    } else if (tmp == "std_dev") {
      trigger_type_ = TriggerMetric::standard_deviation;
    } else if ( tmp == "rel_err") {
      trigger_type_ = TriggerMetric::relative_error;
    } else {
      fatal_error(fmt::format(
        "Invalid volume calculation trigger type '{}' provided.", tmp));
    }

  }

  // Ensure there are no duplicates by copying elements to a set and then
  // comparing the length with the original vector
  std::unordered_set<int> unique_ids(domain_ids_.cbegin(), domain_ids_.cend());
  if (unique_ids.size() != domain_ids_.size()) {
    throw std::runtime_error{"Domain IDs for a volume calculation "
      "must be unique."};
  }

}

std::vector<VolumeCalculation::Result> VolumeCalculation::execute() const
{
  // Shared data that is collected from all threads
  int n = domain_ids_.size();
  std::vector<std::vector<int>> master_indices(n); // List of material indices for each domain
  std::vector<std::vector<int>> master_hits(n); // Number of hits for each material in each domain
  int iterations = 0;

  // Divide work over MPI processes
  size_t min_samples = n_samples_ / mpi::n_procs;
  size_t remainder = n_samples_ % mpi::n_procs;
  size_t i_start, i_end;
  if (mpi::rank < remainder) {
    i_start = (min_samples + 1)*mpi::rank;
    i_end = i_start + min_samples + 1;
  } else {
    i_start = (min_samples + 1)*remainder + (mpi::rank - remainder)*min_samples;
    i_end = i_start + min_samples;
  }

  while (true) {

    #pragma omp parallel
    {
      // Variables that are private to each thread
      std::vector<std::vector<int>> indices(n);
      std::vector<std::vector<int>> hits(n);
      Particle p;

      // Sample locations and count hits
      #pragma omp for
      for (size_t i = i_start; i < i_end; i++) {
        int64_t id = iterations * n_samples_ + i;
        uint64_t seed = init_seed(id, STREAM_VOLUME);

        p.n_coord_ = 1;
        Position xi {prn(&seed), prn(&seed), prn(&seed)};
        p.r() = lower_left_ + xi*(upper_right_ - lower_left_);
        p.u() = {0.5, 0.5, 0.5};

        // If this location is not in the geometry at all, move on to next block
        if (!find_cell(&p, false)) continue;

        if (domain_type_ == TallyDomain::MATERIAL) {
          if (p.material_ != MATERIAL_VOID) {
            for (int i_domain = 0; i_domain < n; i_domain++) {
              if (model::materials[p.material_]->id_ == domain_ids_[i_domain]) {
                this->check_hit(p.material_, indices[i_domain], hits[i_domain]);
                break;
              }
            }
          }
        } else if (domain_type_ == TallyDomain::CELL) {
          for (int level = 0; level < p.n_coord_; ++level) {
            for (int i_domain=0; i_domain < n; i_domain++) {
              if (model::cells[p.coord_[level].cell]->id_ == domain_ids_[i_domain]) {
                this->check_hit(p.material_, indices[i_domain], hits[i_domain]);
                break;
              }
            }
          }
        } else if (domain_type_ == TallyDomain::UNIVERSE) {
          for (int level = 0; level < p.n_coord_; ++level) {
            for (int i_domain = 0; i_domain < n; ++i_domain) {
              if (model::universes[p.coord_[level].universe]->id_ == domain_ids_[i_domain]) {
                check_hit(p.material_, indices[i_domain], hits[i_domain]);
                break;
              }
            }
          }
        }
      }

      // At this point, each thread has its own pair of index/hits lists and we now
      // need to reduce them. OpenMP is not nearly smart enough to do this on its own,
      // so we have to manually reduce them

  #ifdef _OPENMP
      int n_threads = omp_get_num_threads();
  #else
      int n_threads = 1;
  #endif

      #pragma omp for ordered schedule(static)
      for (int i = 0; i < n_threads; ++i) {
        #pragma omp ordered
        for (int i_domain = 0; i_domain < n; ++i_domain) {
          for (int j = 0; j < indices[i_domain].size(); ++j) {
            // Check if this material has been added to the master list and if so,
            // accumulate the number of hits
            bool already_added = false;
            for (int k = 0; k < master_indices[i_domain].size(); k++) {
              if (indices[i_domain][j] == master_indices[i_domain][k]) {
                master_hits[i_domain][k] += hits[i_domain][j];
                already_added = true;
              }
            }
            if (!already_added) {
              // If we made it here, the material hasn't yet been added to the master
              // list, so add entries to the master indices and master hits lists
              master_indices[i_domain].push_back(indices[i_domain][j]);
              master_hits[i_domain].push_back(hits[i_domain][j]);
            }
          }
        }
      }
    } // omp parallel

    // Reduce hits onto master process

    // Determine volume of bounding box
    Position d {upper_right_ - lower_left_};
    double volume_sample = d.x*d.y*d.z;

    // bump iteration counter and get total number
    // of samples at this point
    iterations++;
    size_t total_samples = iterations * n_samples_;

    // reset
    double trigger_val = -INFTY;

    // Set size for members of the Result struct
    std::vector<Result> results(n);

    for (int i_domain = 0; i_domain < n; ++i_domain) {
      // Get reference to result for this domain
      auto& result {results[i_domain]};

      // Create 2D array to store atoms/uncertainty for each nuclide. Later this
      // is compressed into vectors storing only those nuclides that are non-zero
      auto n_nuc = data::nuclides.size();
      xt::xtensor<double, 2> atoms({n_nuc, 2}, 0.0);

#ifdef OPENMC_MPI
      if (mpi::master) {
        for (int j = 1; j < mpi::n_procs; j++) {
          int q;
          MPI_Recv(&q, 1, MPI_INTEGER, j, 0, mpi::intracomm, MPI_STATUS_IGNORE);
          int buffer[2*q];
          MPI_Recv(&buffer[0], 2*q, MPI_INTEGER, j, 1, mpi::intracomm, MPI_STATUS_IGNORE);
          for (int k = 0; k < q; ++k) {
            for (int m = 0; m < master_indices[i_domain].size(); ++m) {
              if (buffer[2*k] == master_indices[i_domain][m]) {
                master_hits[i_domain][m] += buffer[2*k + 1];
                break;
              }
            }
          }
        }
      } else {
        int q = master_indices[i_domain].size();
        int buffer[2*q];
        for (int k = 0; k < q; ++k) {
          buffer[2*k] = master_indices[i_domain][k];
          buffer[2*k + 1] = master_hits[i_domain][k];
        }

        MPI_Send(&q, 1, MPI_INTEGER, 0, 0, mpi::intracomm);
        MPI_Send(&buffer[0], 2*q, MPI_INTEGER, 0, 1, mpi::intracomm);
      }
#endif

      if (mpi::master) {
        int total_hits = 0;
        for (int j = 0; j < master_indices[i_domain].size(); ++j) {
          total_hits += master_hits[i_domain][j];
          double f = static_cast<double>(master_hits[i_domain][j]) / total_samples;
          double var_f = f*(1.0 - f) / total_samples;

          int i_material = master_indices[i_domain][j];
          if (i_material == MATERIAL_VOID) continue;

          const auto& mat = model::materials[i_material];
          for (int k = 0; k < mat->nuclide_.size(); ++k) {
            // Accumulate nuclide density
            int i_nuclide = mat->nuclide_[k];
            atoms(i_nuclide, 0) += mat->atom_density_[k] * f;
            atoms(i_nuclide, 1) += std::pow(mat->atom_density_[k], 2) * var_f;
          }
        }

        // Determine volume
        result.volume[0] = static_cast<double>(total_hits) / total_samples * volume_sample;
        result.volume[1] = std::sqrt(result.volume[0]
          * (volume_sample - result.volume[0]) / total_samples);
        result.iterations = iterations;

        // update threshold value if needed
        if (trigger_type_ != TriggerMetric::not_active) {
          double val = 0.0;
          switch (trigger_type_) {
            case TriggerMetric::standard_deviation:
              val = result.volume[1];
              break;
            case TriggerMetric::relative_error:
              val = result.volume[0] == 0.0 ? INFTY : result.volume[1] / result.volume[0];
              break;
            case TriggerMetric::variance:
              val = result.volume[1] * result.volume[1];
              break;
            default:
              break;
          }
          // update max if entry is valid
          if (val > 0.0) { trigger_val = std::max(trigger_val, val); }
        }

        for (int j = 0; j < n_nuc; ++j) {
          // Determine total number of atoms. At this point, we have values in
          // atoms/b-cm. To get to atoms we multiply by 10^24 V.
          double mean = 1.0e24 * volume_sample * atoms(j, 0);
          double stdev = 1.0e24 * volume_sample * std::sqrt(atoms(j, 1));

          // Convert full arrays to vectors
          if (mean > 0.0) {
            result.nuclides.push_back(j);
            result.atoms.push_back(mean);
            result.uncertainty.push_back(stdev);
          }
        }
      }
    } // end domain loop

    // if no trigger is applied, we're done
    if (trigger_type_ == TriggerMetric::not_active) { return results; }

#ifdef OPENMC_MPI
    // update maximum error value on all processes
    MPI_Bcast(&trigger_val, 1, MPI_DOUBLE, 0, mpi::intracomm);
#endif

    // return results of the calculation
    if (trigger_val < threshold_) { return results; }

#ifdef OPENMC_MPI
    // if iterating in an MPI run, need to zero indices and hits so they aren't counted twice
    if (!mpi::master) {
      for (auto& v : master_indices) { std::fill(v.begin(), v.end(), 0); }
      for (auto& v : master_hits) { std::fill(v.begin(), v.end(), 0); }
    }
#endif

  } // end while
}

void VolumeCalculation::to_hdf5(const std::string& filename,
  const std::vector<Result>& results) const
{
  // Create HDF5 file
  hid_t file_id = file_open(filename, 'w');

  // Write header info
  write_attribute(file_id, "filetype", "volume");
  write_attribute(file_id, "version", VERSION_VOLUME);
  write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  write_attribute(file_id, "date_and_time", time_stamp());

  // Write basic metadata
  write_attribute(file_id, "samples", n_samples_);
  write_attribute(file_id, "lower_left", lower_left_);
  write_attribute(file_id, "upper_right", upper_right_);
  // Write trigger info
  if (trigger_type_ != TriggerMetric::not_active) {
    write_attribute(file_id, "iterations", results[0].iterations);
    write_attribute(file_id, "threshold", threshold_);
    std::string trigger_str;
    switch(trigger_type_) {
      case TriggerMetric::variance:
        trigger_str = "variance";
        break;
      case TriggerMetric::standard_deviation:
        trigger_str = "std_dev";
        break;
      case TriggerMetric::relative_error:
        trigger_str = "rel_err";
        break;
      default:
        break;
    }
    write_attribute(file_id, "trigger_type", trigger_str);
  } else {
    write_attribute(file_id, "iterations", 1);
  }

  if (domain_type_ == TallyDomain::CELL) {
    write_attribute(file_id, "domain_type", "cell");
  }
  else if (domain_type_ == TallyDomain::MATERIAL) {
    write_attribute(file_id, "domain_type", "material");
  }
  else if (domain_type_ == TallyDomain::UNIVERSE) {
    write_attribute(file_id, "domain_type", "universe");
  }

  for (int i = 0; i < domain_ids_.size(); ++i)
  {
    hid_t group_id = create_group(file_id, fmt::format("domain_{}", domain_ids_[i]));

    // Write volume for domain
    const auto& result {results[i]};
    write_dataset(group_id, "volume", result.volume);

    // Create array of nuclide names from the vector
    auto n_nuc = result.nuclides.size();

    std::vector<std::string> nucnames;
    for (int i_nuc : result.nuclides) {
      nucnames.push_back(data::nuclides[i_nuc]->name_);
    }

    // Create array of total # of atoms with uncertainty for each nuclide
    xt::xtensor<double, 2> atom_data({n_nuc, 2});
    xt::view(atom_data, xt::all(), 0) = xt::adapt(result.atoms);
    xt::view(atom_data, xt::all(), 1) = xt::adapt(result.uncertainty);

    // Write results
    write_dataset(group_id, "nuclides", nucnames);
    write_dataset(group_id, "atoms", atom_data);

    close_group(group_id);
  }

  file_close(file_id);
}

void VolumeCalculation::check_hit(int i_material, std::vector<int>& indices,
  std::vector<int>& hits) const
{

  // Check if this material was previously hit and if so, increment count
  bool already_hit = false;
  for (int j = 0; j < indices.size(); j++) {
    if (indices[j] == i_material) {
      hits[j]++;
      already_hit = true;
    }
  }

  // If the material was not previously hit, append an entry to the material
  // indices and hits lists
  if (!already_hit) {
    indices.push_back(i_material);
    hits.push_back(1);
  }
}

void free_memory_volume()
{
  openmc::model::volume_calcs.clear();
}

} // namespace openmc

//==============================================================================
// OPENMC_CALCULATE_VOLUMES runs each of the stochastic volume calculations
// that the user has specified and writes results to HDF5 files
//==============================================================================

int openmc_calculate_volumes() {
  using namespace openmc;

  if (mpi::master) {
    header("STOCHASTIC VOLUME CALCULATION", 3);
  }
  Timer time_volume;
  time_volume.start();

  for (int i = 0; i < model::volume_calcs.size(); ++i) {
    if (mpi::master) {
      write_message(fmt::format("Running volume calculation {}...", i + 1), 4);
    }

    // Run volume calculation
    const auto& vol_calc {model::volume_calcs[i]};
    auto results = vol_calc.execute();

    if (mpi::master) {
      std::string domain_type;
      if (vol_calc.domain_type_ == VolumeCalculation::TallyDomain::CELL) {
        domain_type = "  Cell ";
      } else if (vol_calc.domain_type_ == VolumeCalculation::TallyDomain::MATERIAL) {
        domain_type = "  Material ";
      } else {
        domain_type = "  Universe ";
      }

      // Display domain volumes
      for (int j = 0; j < vol_calc.domain_ids_.size(); j++) {
        write_message(fmt::format("{}{}: {} +/- {} cm^3", domain_type,
          vol_calc.domain_ids_[j], results[j].volume[0], results[j].volume[1]), 4);
      }

      // Write volumes to HDF5 file
      std::string filename = fmt::format("{}volume_{}.h5",
        settings::path_output, i + 1);
      vol_calc.to_hdf5(filename, results);
    }

  }

  // Show elapsed time
  time_volume.stop();
  if (mpi::master) {
    write_message(fmt::format("Elapsed time: {} s", time_volume.elapsed()), 6);
  }

  return 0;
}
