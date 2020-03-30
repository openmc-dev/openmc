#include "openmc/state_point.h"

#include <algorithm>
#include <cstdint> // for int64_t
#include <string>
#include <vector>

#include <fmt/core.h>
#include "xtensor/xbuilder.hpp" // for empty_like
#include "xtensor/xview.hpp"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/tally.h"
#include "openmc/timer.h"

namespace openmc {

extern "C" int
openmc_statepoint_write(const char* filename, bool* write_source)
{
  // Set the filename
  std::string filename_;
  if (filename) {
    filename_ = filename;
  } else {
    // Determine width for zero padding
    int w = std::to_string(settings::n_max_batches).size();

    // Set filename for state point
    filename_ = fmt::format("{0}statepoint.{1:0{2}}.h5",
      settings::path_output, simulation::current_batch, w);
  }

  // Determine whether or not to write the source bank
  bool write_source_ = write_source ? *write_source : true;

  // Write message
  write_message("Creating state point " + filename_ + "...", 5);

  hid_t file_id;
  if (mpi::master) {
    // Create statepoint file
    file_id = file_open(filename_, 'w');

    // Write file type
    write_attribute(file_id, "filetype", "statepoint");

    // Write revision number for state point file
    write_attribute(file_id, "version", VERSION_STATEPOINT);

    // Write OpenMC version
    write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
    write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

    // Write current date and time
    write_attribute(file_id, "date_and_time", time_stamp());

    // Write path to input
    write_attribute(file_id, "path", settings::path_input);

    // Write out random number seed
    write_dataset(file_id, "seed", openmc_get_seed());

    // Write run information
    write_dataset(file_id, "energy_mode", settings::run_CE ?
      "continuous-energy" : "multi-group");
    switch (settings::run_mode) {
    case RunMode::FIXED_SOURCE:
      write_dataset(file_id, "run_mode", "fixed source");
      break;
    case RunMode::EIGENVALUE:
      write_dataset(file_id, "run_mode", "eigenvalue");
      break;
    default:
      break;
    }
    write_attribute(file_id, "photon_transport", settings::photon_transport);
    write_dataset(file_id, "n_particles", settings::n_particles);
    write_dataset(file_id, "n_batches", settings::n_batches);

    // Write out current batch number
    write_dataset(file_id, "current_batch", simulation::current_batch);

    // Indicate whether source bank is stored in statepoint
    write_attribute(file_id, "source_present", write_source_);

    // Write out information for eigenvalue run
    if (settings::run_mode == RunMode::EIGENVALUE)
      write_eigenvalue_hdf5(file_id);

    hid_t tallies_group = create_group(file_id, "tallies");

    // Write meshes
    meshes_to_hdf5(tallies_group);

    // Write information for derivatives
    if (!model::tally_derivs.empty()) {
      hid_t derivs_group = create_group(tallies_group, "derivatives");
      for (const auto& deriv : model::tally_derivs) {
        hid_t deriv_group = create_group(derivs_group,
          "derivative " + std::to_string(deriv.id));
        write_dataset(deriv_group, "material", deriv.diff_material);
        if (deriv.variable == DerivativeVariable::DENSITY) {
          write_dataset(deriv_group, "independent variable", "density");
        } else if (deriv.variable == DerivativeVariable::NUCLIDE_DENSITY) {
          write_dataset(deriv_group, "independent variable", "nuclide_density");
          write_dataset(deriv_group, "nuclide",
            data::nuclides[deriv.diff_nuclide]->name_);
        } else if (deriv.variable == DerivativeVariable::TEMPERATURE) {
          write_dataset(deriv_group, "independent variable", "temperature");
        } else {
          fatal_error("Independent variable for derivative "
            + std::to_string(deriv.id) + " not defined in state_point.cpp");
        }
        close_group(deriv_group);
      }
      close_group(derivs_group);
    }

    // Write information for filters
    hid_t filters_group = create_group(tallies_group, "filters");
    write_attribute(filters_group, "n_filters", model::tally_filters.size());
    if (!model::tally_filters.empty()) {
      // Write filter IDs
      std::vector<int32_t> filter_ids;
      filter_ids.reserve(model::tally_filters.size());
      for (const auto& filt : model::tally_filters)
        filter_ids.push_back(filt->id());
      write_attribute(filters_group, "ids", filter_ids);

      // Write info for each filter
      for (const auto& filt : model::tally_filters) {
        hid_t filter_group = create_group(filters_group,
          "filter " + std::to_string(filt->id()));
        filt->to_statepoint(filter_group);
        close_group(filter_group);
      }
    }
    close_group(filters_group);

    // Write information for tallies
    write_attribute(tallies_group, "n_tallies", model::tallies.size());
    if (!model::tallies.empty()) {
      // Write tally IDs
      std::vector<int32_t> tally_ids;
      tally_ids.reserve(model::tallies.size());
      for (const auto& tally : model::tallies)
        tally_ids.push_back(tally->id_);
      write_attribute(tallies_group, "ids", tally_ids);

#ifdef DAGMC
      // write unstructured mesh tallies to VTK if possible
      write_unstructured_mesh_results();
#endif

      // Write all tally information except results
      for (const auto& tally : model::tallies) {
        hid_t tally_group = create_group(tallies_group,
          "tally " + std::to_string(tally->id_));

        write_dataset(tally_group, "name",  tally->name_);

        if (tally->writable_) {
          write_attribute(tally_group, "internal", 0);
        } else {
          write_attribute(tally_group, "internal", 1);
          close_group(tally_group);
          continue;
        }

        if (tally->estimator_ == TallyEstimator::ANALOG) {
          write_dataset(tally_group, "estimator", "analog");
        } else if (tally->estimator_ == TallyEstimator::TRACKLENGTH) {
          write_dataset(tally_group, "estimator", "tracklength");
        } else if (tally->estimator_ == TallyEstimator::COLLISION) {
          write_dataset(tally_group, "estimator", "collision");
        }

        write_dataset(tally_group, "n_realizations", tally->n_realizations_);

        // Write the ID of each filter attached to this tally
        write_dataset(tally_group, "n_filters", tally->filters().size());
        if (!tally->filters().empty()) {
          std::vector<int32_t> filter_ids;
          filter_ids.reserve(tally->filters().size());
          for (auto i_filt : tally->filters())
            filter_ids.push_back(model::tally_filters[i_filt]->id());
          write_dataset(tally_group, "filters", filter_ids);
        }

        // Write the nuclides this tally scores
        std::vector<std::string> nuclides;
        for (auto i_nuclide : tally->nuclides_) {
          if (i_nuclide == -1) {
            nuclides.push_back("total");
          } else {
            if (settings::run_CE) {
              nuclides.push_back(data::nuclides[i_nuclide]->name_);
            } else {
              nuclides.push_back(data::mg.nuclides_[i_nuclide].name);
            }
          }
        }
        write_dataset(tally_group, "nuclides", nuclides);

        if (tally->deriv_ != C_NONE) write_dataset(tally_group, "derivative",
          model::tally_derivs[tally->deriv_].id);

        // Write the tally score bins
        std::vector<std::string> scores;
        for (auto sc : tally->scores_) scores.push_back(reaction_name(sc));
        write_dataset(tally_group, "n_score_bins", scores.size());
        write_dataset(tally_group, "score_bins", scores);

        close_group(tally_group);
      }

    }

    if (settings::reduce_tallies) {
      // Write global tallies
      write_dataset(file_id, "global_tallies", simulation::global_tallies);

      // Write tallies
      if (model::active_tallies.size() > 0) {
        // Indicate that tallies are on
        write_attribute(file_id, "tallies_present", 1);

        // Write all tally results
        for (const auto& tally : model::tallies) {
          if (!tally->writable_) continue;
          // Write sum and sum_sq for each bin
          std::string name = "tally " + std::to_string(tally->id_);
          hid_t tally_group = open_group(tallies_group, name.c_str());
          auto& results = tally->results_;
          write_tally_results(tally_group, results.shape()[0],
            results.shape()[1], results.data());
          close_group(tally_group);
        }
      } else {
        // Indicate tallies are off
        write_attribute(file_id, "tallies_present", 0);
      }
    }

    close_group(tallies_group);
  }

  // Check for the no-tally-reduction method
  if (!settings::reduce_tallies) {
    // If using the no-tally-reduction method, we need to collect tally
    // results before writing them to the state point file.
    write_tally_results_nr(file_id);

  } else if (mpi::master) {
    // Write number of global realizations
    write_dataset(file_id, "n_realizations", simulation::n_realizations);

    // Write out the runtime metrics.
    using namespace simulation;
    hid_t runtime_group = create_group(file_id, "runtime");
    write_dataset(runtime_group, "total initialization",  time_initialize.elapsed());
    write_dataset(runtime_group, "reading cross sections", time_read_xs.elapsed());
    write_dataset(runtime_group, "simulation", time_inactive.elapsed()
      + time_active.elapsed());
    write_dataset(runtime_group, "transport", time_transport.elapsed());
    if (settings::run_mode == RunMode::EIGENVALUE) {
      write_dataset(runtime_group, "inactive batches", time_inactive.elapsed());
    }
    write_dataset(runtime_group, "active batches", time_active.elapsed());
    if (settings::run_mode == RunMode::EIGENVALUE) {
      write_dataset(runtime_group, "synchronizing fission bank", time_bank.elapsed());
      write_dataset(runtime_group, "sampling source sites", time_bank_sample.elapsed());
      write_dataset(runtime_group, "SEND-RECV source sites", time_bank_sendrecv.elapsed());
    } else {
      write_dataset(runtime_group, "sampling source sites", time_sample_source.elapsed());
    }
    write_dataset(runtime_group, "accumulating tallies", time_tallies.elapsed());
    write_dataset(runtime_group, "total", time_total.elapsed());
    close_group(runtime_group);

    file_close(file_id);
  }

#ifdef PHDF5
  bool parallel = true;
#else
  bool parallel = false;
#endif

  // Write the source bank if desired
  if (write_source_) {
    if (mpi::master || parallel) file_id = file_open(filename_, 'a', true);
    write_source_bank(file_id);
    if (mpi::master || parallel) file_close(file_id);
  }

  return 0;
}

void restart_set_keff()
{
  if (simulation::restart_batch > settings::n_inactive) {
    for (int i = settings::n_inactive; i < simulation::restart_batch; ++i) {
      simulation::k_sum[0] += simulation::k_generation[i];
      simulation::k_sum[1] += std::pow(simulation::k_generation[i], 2);
    }
    int n = settings::gen_per_batch*simulation::n_realizations;
    simulation::keff = simulation::k_sum[0] / n;
  } else {
    simulation::keff = simulation::k_generation.back();
  }
}

void load_state_point()
{
  // Write message
  write_message("Loading state point " + settings::path_statepoint + "...", 5);

  // Open file for reading
  hid_t file_id = file_open(settings::path_statepoint.c_str(), 'r', true);

  // Read filetype
  std::string word;
  read_attribute(file_id, "filetype", word);
  if (word != "statepoint") {
    fatal_error("OpenMC tried to restart from a non-statepoint file.");
  }

  // Read revision number for state point file and make sure it matches with
  // current version
  std::array<int, 2> array;
  read_attribute(file_id, "version", array);
  if (array != VERSION_STATEPOINT) {
    fatal_error("State point version does not match current version in OpenMC.");
  }

  // Read and overwrite random number seed
  int64_t seed;
  read_dataset(file_id, "seed", seed);
  openmc_set_seed(seed);

  // It is not impossible for a state point to be generated from a CE run but
  // to be loaded in to an MG run (or vice versa), check to prevent that.
  read_dataset(file_id, "energy_mode", word);
  if (word == "multi-group" && settings::run_CE) {
    fatal_error("State point file is from multigroup run but current run is "
      "continous energy.");
  } else if (word == "continuous-energy" && !settings::run_CE) {
    fatal_error("State point file is from continuous-energy run but current "
      "run is multigroup!");
  }

  // Read and overwrite run information except number of batches
  read_dataset(file_id, "run_mode", word);
  if (word == "fixed source") {
    settings::run_mode = RunMode::FIXED_SOURCE;
  } else if (word == "eigenvalue") {
    settings::run_mode = RunMode::EIGENVALUE;
  }
  read_attribute(file_id, "photon_transport", settings::photon_transport);
  read_dataset(file_id, "n_particles", settings::n_particles);
  int temp;
  read_dataset(file_id, "n_batches", temp);

  // Take maximum of statepoint n_batches and input n_batches
  settings::n_batches = std::max(settings::n_batches, temp);

  // Read batch number to restart at
  read_dataset(file_id, "current_batch", simulation::restart_batch);

  if (simulation::restart_batch > settings::n_batches) {
    fatal_error("The number batches specified in settings.xml is fewer "
      " than the number of batches in the given statepoint file.");
  }

  // Logical flag for source present in statepoint file
  bool source_present;
  read_attribute(file_id, "source_present", source_present);

  // Read information specific to eigenvalue run
  if (settings::run_mode == RunMode::EIGENVALUE) {
    read_dataset(file_id, "n_inactive", temp);
    read_eigenvalue_hdf5(file_id);

    // Take maximum of statepoint n_inactive and input n_inactive
    settings::n_inactive = std::max(settings::n_inactive, temp);

    // Check to make sure source bank is present
    if (settings::path_sourcepoint == settings::path_statepoint &&
        !source_present) {
      fatal_error("Source bank must be contained in statepoint restart file");
    }
  }

  // Read number of realizations for global tallies
  read_dataset(file_id, "n_realizations", simulation::n_realizations);

  // Set k_sum, keff, and current_batch based on whether restart file is part
  // of active cycle or inactive cycle
  if (settings::run_mode == RunMode::EIGENVALUE) {
    restart_set_keff();
  }

  // Set current batch number
  simulation::current_batch = simulation::restart_batch;

  // Read tallies to master. If we are using Parallel HDF5, all processes
  // need to be included in the HDF5 calls.
#ifdef PHDF5
  if (true) {
#else
  if (mpi::master) {
#endif
    // Read global tally data
    read_dataset_lowlevel(file_id, "global_tallies", H5T_NATIVE_DOUBLE,
      H5S_ALL, false, simulation::global_tallies.data());

    // Check if tally results are present
    bool present;
    read_attribute(file_id, "tallies_present", present);

    // Read in sum and sum squared
    if (present) {
      hid_t tallies_group = open_group(file_id, "tallies");

      for (auto& tally : model::tallies) {
        // Read sum, sum_sq, and N for each bin
        std::string name = "tally " + std::to_string(tally->id_);
        hid_t tally_group = open_group(tallies_group, name.c_str());

        int internal=0;
        if (attribute_exists(tally_group, "internal")) {
          read_attribute(tally_group, "internal", internal);
        }
        if (internal) {
          tally->writable_ = false;
        } else {

          auto& results = tally->results_;
          read_tally_results(tally_group, results.shape()[0],
            results.shape()[1], results.data());
          read_dataset(tally_group, "n_realizations", tally->n_realizations_);
          close_group(tally_group);
        }
      }

      close_group(tallies_group);
    }
  }

  // Read source if in eigenvalue mode
  if (settings::run_mode == RunMode::EIGENVALUE) {

    // Check if source was written out separately
    if (!source_present) {

      // Close statepoint file
      file_close(file_id);

      // Write message
      write_message("Loading source file " + settings::path_sourcepoint
        + "...", 5);

      // Open source file
      file_id = file_open(settings::path_source.c_str(), 'r', true);
    }

    // Read source
    read_source_bank(file_id);

  }

  // Close file
  file_close(file_id);
}


hid_t h5banktype() {
  // Create compound type for position
  hid_t postype = H5Tcreate(H5T_COMPOUND, sizeof(struct Position));
  H5Tinsert(postype, "x", HOFFSET(Position, x), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "y", HOFFSET(Position, y), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "z", HOFFSET(Position, z), H5T_NATIVE_DOUBLE);

  // Create bank datatype
  hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(struct Particle::Bank));
  H5Tinsert(banktype, "r", HOFFSET(Particle::Bank, r), postype);
  H5Tinsert(banktype, "u", HOFFSET(Particle::Bank, u), postype);
  H5Tinsert(banktype, "E", HOFFSET(Particle::Bank, E), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "wgt", HOFFSET(Particle::Bank, wgt), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "delayed_group", HOFFSET(Particle::Bank, delayed_group), H5T_NATIVE_INT);
  H5Tinsert(banktype, "particle", HOFFSET(Particle::Bank, particle), H5T_NATIVE_INT);

  H5Tclose(postype);
  return banktype;
}

void
write_source_point(const char* filename)
{
  // When using parallel HDF5, the file is written to collectively by all
  // processes. With MPI-only, the file is opened and written by the master
  // (note that the call to write_source_bank is by all processes since slave
  // processes need to send source bank data to the master.
#ifdef PHDF5
  bool parallel = true;
#else
  bool parallel = false;
#endif

  std::string filename_;
  if (filename) {
    filename_ = filename;
  } else {
    // Determine width for zero padding
    int w = std::to_string(settings::n_max_batches).size();

    filename_ = fmt::format("{0}source.{1:0{2}}.h5",
      settings::path_output, simulation::current_batch, w);
  }

  hid_t file_id;
  if (mpi::master || parallel) {
    file_id = file_open(filename_, 'w', true);
    write_attribute(file_id, "filetype", "source");
  }

  // Get pointer to source bank and write to file
  write_source_bank(file_id);

  if (mpi::master || parallel) file_close(file_id);
}

void
write_source_bank(hid_t group_id)
{
  hid_t banktype = h5banktype();

#ifdef PHDF5
  // Set size of total dataspace for all procs and rank
  hsize_t dims[] {static_cast<hsize_t>(settings::n_particles)};
  hid_t dspace = H5Screate_simple(1, dims, nullptr);
  hid_t dset = H5Dcreate(group_id, "source_bank", banktype, dspace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create another data space but for each proc individually
  hsize_t count[] {static_cast<hsize_t>(simulation::work_per_rank)};
  hid_t memspace = H5Screate_simple(1, count, nullptr);

  // Select hyperslab for this dataspace
  hsize_t start[] {static_cast<hsize_t>(simulation::work_index[mpi::rank])};
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Set up the property list for parallel writing
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

  // Write data to file in parallel
  H5Dwrite(dset, banktype, memspace, dspace, plist, simulation::source_bank.data());

  // Free resources
  H5Sclose(dspace);
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Pclose(plist);

#else

  if (mpi::master) {
    // Create dataset big enough to hold all source sites
    hsize_t dims[] {static_cast<hsize_t>(settings::n_particles)};
    hid_t dspace = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate(group_id, "source_bank", banktype, dspace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Save source bank sites since the souce_bank array is overwritten below
#ifdef OPENMC_MPI
    std::vector<Particle::Bank> temp_source {simulation::source_bank.begin(),
      simulation::source_bank.begin() + simulation::work_per_rank};
#endif

    for (int i = 0; i < mpi::n_procs; ++i) {
      // Create memory space
      hsize_t count[] {static_cast<hsize_t>(simulation::work_index[i+1] -
        simulation::work_index[i])};
      hid_t memspace = H5Screate_simple(1, count, nullptr);

#ifdef OPENMC_MPI
      // Receive source sites from other processes
      if (i > 0)
        MPI_Recv(simulation::source_bank.data(), count[0], mpi::bank, i, i,
                 mpi::intracomm, MPI_STATUS_IGNORE);
#endif

      // Select hyperslab for this dataspace
      dspace = H5Dget_space(dset);
      hsize_t start[] {static_cast<hsize_t>(simulation::work_index[i])};
      H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Write data to hyperslab
      H5Dwrite(dset, banktype, memspace, dspace, H5P_DEFAULT,
        simulation::source_bank.data());

      H5Sclose(memspace);
      H5Sclose(dspace);
    }

    // Close all ids
    H5Dclose(dset);

#ifdef OPENMC_MPI
    // Restore state of source bank
    std::copy(temp_source.begin(), temp_source.end(), simulation::source_bank.begin());
#endif
  } else {
#ifdef OPENMC_MPI
    MPI_Send(simulation::source_bank.data(), simulation::work_per_rank, mpi::bank,
      0, mpi::rank, mpi::intracomm);
#endif
  }
#endif

  H5Tclose(banktype);
}


void read_source_bank(hid_t group_id)
{
  hid_t banktype = h5banktype();

  // Open the dataset
  hid_t dset = H5Dopen(group_id, "source_bank", H5P_DEFAULT);

  // Create another data space but for each proc individually
  hsize_t dims[] {static_cast<hsize_t>(simulation::work_per_rank)};
  hid_t memspace = H5Screate_simple(1, dims, nullptr);

  // Make sure source bank is big enough
  hid_t dspace = H5Dget_space(dset);
  hsize_t dims_all[1];
  H5Sget_simple_extent_dims(dspace, dims_all, nullptr);
  if (simulation::work_index[mpi::n_procs] > dims_all[0]) {
    fatal_error("Number of source sites in source file is less "
                "than number of source particles per generation.");
  }

  // Select hyperslab for each process
  hsize_t start[] {static_cast<hsize_t>(simulation::work_index[mpi::rank])};
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, dims, nullptr);

#ifdef PHDF5
    // Read data in parallel
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset, banktype, memspace, dspace, plist, simulation::source_bank.data());
  H5Pclose(plist);
#else
  H5Dread(dset, banktype, memspace, dspace, H5P_DEFAULT, simulation::source_bank.data());
#endif

  // Close all ids
  H5Sclose(dspace);
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Tclose(banktype);
}

#ifdef DAGMC
void write_unstructured_mesh_results() {

  for (auto& tally : model::tallies) {
    for (auto filter_idx : tally->filters()) {
      auto& filter = model::tally_filters[filter_idx];
      if (filter->type() != "mesh") continue;

      // check if the filter uses an unstructured mesh
      auto mesh_filter = dynamic_cast<MeshFilter*>(filter.get());
      auto mesh_idx = mesh_filter->mesh();
      auto umesh = dynamic_cast<UnstructuredMesh*>(model::meshes[mesh_idx].get());

      if (!umesh) continue;

      // if this tally has more than one filter, print
      // warning and skip writing the mesh
      if (tally->filters().size() > 1) {
        warning(fmt::format("Skipping unstructured mesh writing for tally "
                            "{}. More than one filter is present on the tally.",
                            tally->id_));
        break;
      }

      int n_realizations = tally->n_realizations_;

      // write each score/nuclide combination for this tally
      for (int i_score = 0; i_score < tally->scores_.size(); i_score++) {
        for (int i_nuc = 0; i_nuc < tally->nuclides_.size(); i_nuc++) {

          // index for this nuclide and score
          int nuc_score_idx = i_score + i_nuc*tally->scores_.size();

          // construct result vectors
          std::vector<double> mean_vec, std_dev_vec;
          for (int j = 0; j < tally->results_.shape()[0]; j++) {
            // mean
            double mean = tally->results_(j, nuc_score_idx, TallyResult::SUM) / n_realizations;
            mean_vec.push_back(mean);
            // std. dev.
            double sum_sq = tally->results_(j , nuc_score_idx, TallyResult::SUM_SQ);
            if (n_realizations > 1) {
              double std_dev = sum_sq/n_realizations - mean*mean;
              std_dev = std::sqrt(std_dev / (n_realizations - 1));
              std_dev_vec.push_back(std_dev);
            } else {
              std_dev_vec.push_back(0.0);
            }
          }

          // generate a name for the value
          std::string nuclide_name = "total"; // start with total by default
          if (tally->nuclides_[i_nuc] > -1) {
            nuclide_name = data::nuclides[tally->nuclides_[i_nuc]]->name_;
          }

          std::string score_name = tally->score_name(i_score);
          auto score_str = fmt::format("{}_{}", score_name, nuclide_name);
          umesh->set_score_data(score_str, mean_vec, std_dev_vec);
        }
      }

      // Generate a file name based on the tally id
      // and the current batch number
      int w = std::to_string(settings::n_max_batches).size();
      std::string filename = fmt::format("tally_{0}.{1:0{2}}",
                                         tally->id_,
                                         simulation::current_batch,
                                         w);
      // Write the unstructured mesh and data to file
      umesh->write(filename);
    }
  }
}
#endif

void write_tally_results_nr(hid_t file_id)
{
  // ==========================================================================
  // COLLECT AND WRITE GLOBAL TALLIES

  hid_t tallies_group;
  if (mpi::master) {
    // Write number of realizations
    write_dataset(file_id, "n_realizations", simulation::n_realizations);

    // Write number of global tallies
    write_dataset(file_id, "n_global_tallies", N_GLOBAL_TALLIES);

    tallies_group = open_group(file_id, "tallies");
  }

  // Get global tallies
  auto& gt = simulation::global_tallies;

#ifdef OPENMC_MPI
  // Reduce global tallies
  xt::xtensor<double, 2> gt_reduced = xt::empty_like(gt);
  MPI_Reduce(gt.data(), gt_reduced.data(), gt.size(), MPI_DOUBLE,
    MPI_SUM, 0, mpi::intracomm);

  // Transfer values to value on master
  if (mpi::master) {
    if (simulation::current_batch == settings::n_max_batches ||
        simulation::satisfy_triggers) {
      std::copy(gt_reduced.begin(), gt_reduced.end(), gt.begin());
    }
  }
#endif

  // Write out global tallies sum and sum_sq
  if (mpi::master) {
    write_dataset(file_id, "global_tallies", gt);
  }

  for (const auto& t : model::tallies) {
    // Skip any tallies that are not active
    if (!t->active_) continue;
    if (!t->writable_) continue;

    if (mpi::master && !object_exists(file_id, "tallies_present")) {
      write_attribute(file_id, "tallies_present", 1);
    }

    // Get view of accumulated tally values
    auto values_view = xt::view(t->results_, xt::all(), xt::all(),
      xt::range(static_cast<int>(TallyResult::SUM), static_cast<int>(TallyResult::SUM_SQ) + 1));

    // Make copy of tally values in contiguous array
    xt::xtensor<double, 2> values = values_view;

    if (mpi::master) {
      // Open group for tally
      std::string groupname {"tally " + std::to_string(t->id_)};
      hid_t tally_group = open_group(tallies_group, groupname.c_str());

      // The MPI_IN_PLACE specifier allows the master to copy values into
      // a receive buffer without having a temporary variable
#ifdef OPENMC_MPI
      MPI_Reduce(MPI_IN_PLACE, values.data(), values.size(), MPI_DOUBLE,
        MPI_SUM, 0, mpi::intracomm);
#endif

      // At the end of the simulation, store the results back in the
      // regular TallyResults array
      if (simulation::current_batch == settings::n_max_batches ||
          simulation::satisfy_triggers) {
        values_view = values;
      }

      // Put in temporary tally result
      xt::xtensor<double, 3> results_copy = xt::zeros_like(t->results_);
      auto copy_view = xt::view(results_copy, xt::all(), xt::all(),
        xt::range(static_cast<int>(TallyResult::SUM), static_cast<int>(TallyResult::SUM_SQ) + 1));
      copy_view = values;

      // Write reduced tally results to file
      auto shape = results_copy.shape();
      write_tally_results(tally_group, shape[0], shape[1], results_copy.data());

      close_group(tally_group);
    } else {
      // Receive buffer not significant at other processors
#ifdef OPENMC_MPI
      MPI_Reduce(values.data(), nullptr, values.size(), MPI_DOUBLE, MPI_SUM,
            0, mpi::intracomm);
#endif
    }
  }

  if (mpi::master) {
    if (!object_exists(file_id, "tallies_present")) {
      // Indicate that tallies are off
      write_dataset(file_id, "tallies_present", 0);
    }
  }
}

} // namespace openmc
