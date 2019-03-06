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

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"

#include <algorithm> // for copy
#include <cmath> // for pow, sqrt
#include <sstream>
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
    domain_type_ = FILTER_CELL;
  } else if (domain_type == "material") {
    domain_type_ = FILTER_MATERIAL;
  } else if (domain_type == "universe") {
    domain_type_ = FILTER_UNIVERSE;
  } else {
    fatal_error(std::string("Unrecognized domain type for stochastic "
      "volume calculation: " + domain_type));
  }

  // Read domain IDs, bounding corodinates and number of samples
  domain_ids_ = get_node_array<int>(node, "domain_ids");
  lower_left_ = get_node_array<double>(node, "lower_left");
  upper_right_ = get_node_array<double>(node, "upper_right");
  n_samples_ = std::stoi(get_node_value(node, "samples"));

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

  // Divide work over MPI processes
  int min_samples = n_samples_ / mpi::n_procs;
  int remainder = n_samples_ % mpi::n_procs;
  int i_start, i_end;
  if (mpi::rank < remainder) {
    i_start = (min_samples + 1)*mpi::rank;
    i_end = i_start + min_samples + 1;
  } else {
    i_start = (min_samples + 1)*remainder + (mpi::rank - remainder)*min_samples;
    i_end = i_start + min_samples;
  }

  #pragma omp parallel
  {
    // Variables that are private to each thread
    std::vector<std::vector<int>> indices(n);
    std::vector<std::vector<int>> hits(n);
    Particle p;

    prn_set_stream(STREAM_VOLUME);

    // Sample locations and count hits
    #pragma omp for
    for (int i = i_start; i < i_end; i++) {
      set_particle_seed(i);

      p.n_coord_ = 1;
      Position xi {prn(), prn(), prn()};
      p.r() = lower_left_ + xi*(upper_right_ - lower_left_);
      p.u() = {0.5, 0.5, 0.5};

      // If this location is not in the geometry at all, move on to next block
      if (!find_cell(&p, false)) continue;

      if (domain_type_ == FILTER_MATERIAL) {
        if (p.material_ != MATERIAL_VOID) {
          for (int i_domain = 0; i_domain < n; i_domain++) {
            if (model::materials[p.material_]->id_ == domain_ids_[i_domain]) {
              this->check_hit(p.material_, indices[i_domain], hits[i_domain]);
              break;
            }
          }
        }
      } else if (domain_type_ == FILTER_CELL) {
        for (int level = 0; level < p.n_coord_; ++level) {
          for (int i_domain=0; i_domain < n; i_domain++) {
            if (model::cells[p.coord_[level].cell]->id_ == domain_ids_[i_domain]) {
              this->check_hit(p.material_, indices[i_domain], hits[i_domain]);
              break;
            }
          }
        }
      } else if (domain_type_ == FILTER_UNIVERSE) {
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
    #pragma omp for ordered schedule(static)
    for (int i = 0; i < omp_get_num_threads(); ++i) {
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
#else
    master_indices = indices;
    master_hits = hits;
#endif

    prn_set_stream(STREAM_TRACKING);
  } // omp parallel

  // Reduce hits onto master process

  // Determine volume of bounding box
  Position d {upper_right_ - lower_left_};
  double volume_sample = d.x*d.y*d.z;

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
        double f = static_cast<double>(master_hits[i_domain][j]) / n_samples_;
        double var_f = f*(1.0 - f)/n_samples_;

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
      result.volume[0] = static_cast<double>(total_hits) / n_samples_ * volume_sample;
      result.volume[1] = std::sqrt(result.volume[0]
        * (volume_sample - result.volume[0]) / n_samples_);

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
  }

  return results;
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
  if (domain_type_ == FILTER_CELL) {
    write_attribute(file_id, "domain_type", "cell");
  }
  else if (domain_type_ == FILTER_MATERIAL) {
    write_attribute(file_id, "domain_type", "material");
  }
  else if (domain_type_ == FILTER_UNIVERSE) {
    write_attribute(file_id, "domain_type", "universe");
  }

  for (int i = 0; i < domain_ids_.size(); ++i)
  {
    hid_t group_id = create_group(file_id, "domain_"
      + std::to_string(domain_ids_[i]));

    // Write volume for domain
    const auto& result {results[i]};
    write_dataset(group_id, "volume", result.volume);

    // Create array of nuclide names from the vector
    auto n_nuc = result.nuclides.size();

    if (!result.nuclides.empty()) {
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
    }

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
      write_message("Running volume calculation " + std::to_string(i+1) + "...", 4);
    }

    // Run volume calculation
    const auto& vol_calc {model::volume_calcs[i]};
    auto results = vol_calc.execute();

    if (mpi::master) {
      std::string domain_type;
      if (vol_calc.domain_type_ == FILTER_CELL) {
        domain_type = "  Cell ";
      } else if (vol_calc.domain_type_ == FILTER_MATERIAL) {
        domain_type = "  Material ";
      } else {
        domain_type = "  Universe ";
      }

      // Display domain volumes
      for (int j = 0; j < vol_calc.domain_ids_.size(); j++) {
        std::stringstream msg;
        msg << domain_type << vol_calc.domain_ids_[j] << ": " <<
          results[j].volume[0] << " +/- " << results[j].volume[1] << " cm^3";
        write_message(msg, 4);
      }

      // Write volumes to HDF5 file
      std::string filename = settings::path_output + "volume_"
        + std::to_string(i+1) + ".h5";
      vol_calc.to_hdf5(filename, results);
    }

  }

  // Show elapsed time
  time_volume.stop();
  if (mpi::master) {
    write_message("Elapsed time: " + std::to_string(time_volume.elapsed())
      + " s", 6);
  }

  return 0;
}
