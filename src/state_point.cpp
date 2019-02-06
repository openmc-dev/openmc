#include "openmc/state_point.h"

#include <algorithm>
#include <cstdint> // for int64_t
#include <iomanip> // for setfill, setw
#include <string>
#include <vector>

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
#include "openmc/output.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/timer.h"

namespace openmc {

extern "C" void statepoint_write_f(hid_t file_id);
extern "C" void load_state_point_f(hid_t file_id);

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
    std::stringstream ss;
    ss << settings::path_output << "statepoint." << std::setfill('0')
      << std::setw(w) << simulation::current_batch << ".h5";
    filename_ = ss.str();
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
    case RUN_MODE_FIXEDSOURCE:
      write_dataset(file_id, "run_mode", "fixed source");
      break;
    case RUN_MODE_EIGENVALUE:
      write_dataset(file_id, "run_mode", "eigenvalue");
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
    if (settings::run_mode == RUN_MODE_EIGENVALUE) write_eigenvalue_hdf5(file_id);

    // Write tally-related data
    hid_t tallies_group = create_group(file_id, "tallies");
    meshes_to_hdf5(tallies_group);
    statepoint_write_f(file_id);
    close_group(tallies_group);
  }

  // Check for the no-tally-reduction method
  if (!settings::reduce_tallies) {
    // If using the no-tally-reduction method, we need to collect tally
    // results before writing them to the state point file.

    write_tally_results_nr(file_id);

  } else if (mpi::master) {

    // Write number of global realizations
    write_dataset(file_id, "n_realizations", n_realizations);

    // Write out the runtime metrics.
    using namespace simulation;
    hid_t runtime_group = create_group(file_id, "runtime");
    write_dataset(runtime_group, "total initialization",  time_initialize.elapsed());
    write_dataset(runtime_group, "reading cross sections", time_read_xs.elapsed());
    write_dataset(runtime_group, "simulation", time_inactive.elapsed()
      + time_active.elapsed());
    write_dataset(runtime_group, "transport", time_transport.elapsed());
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      write_dataset(runtime_group, "inactive batches", time_inactive.elapsed());
    }
    write_dataset(runtime_group, "active batches", time_active.elapsed());
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      write_dataset(runtime_group, "synchronizing fission bank", time_bank.elapsed());
      write_dataset(runtime_group, "sampling source sites", time_bank_sample.elapsed());
      write_dataset(runtime_group, "SEND-RECV source sites", time_bank_sendrecv.elapsed());
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
    int n = settings::gen_per_batch*n_realizations;
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
    settings::run_mode = RUN_MODE_FIXEDSOURCE;
  } else if (word == "eigenvalue") {
    settings::run_mode = RUN_MODE_EIGENVALUE;
  }
  read_attribute(file_id, "photon_transport", settings::photon_transport);
  read_dataset(file_id, "n_particles", settings::n_particles);
  int temp;
  read_dataset(file_id, "n_batches", temp);

  // Take maximum of statepoint n_batches and input n_batches
  settings::n_batches = std::max(settings::n_batches, temp);

  // Read batch number to restart at
  read_dataset(file_id, "current_batch", simulation::restart_batch);

  // Check for source in statepoint if needed
  bool source_present;
  read_attribute(file_id, "source_present", source_present);

  if (simulation::restart_batch > settings::n_batches) {
    fatal_error("The number batches specified in settings.xml is fewer "
      " than the number of batches in the given statepoint file.");
  }

  // Read information specific to eigenvalue run
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {
    read_dataset(file_id, "n_inactive", temp);
    read_eigenvalue_hdf5(file_id);

    // Take maximum of statepoint n_inactive and input n_inactive
    settings::n_inactive = std::max(settings::n_inactive, temp);
  }

  // Read number of realizations for global tallies
  read_dataset(file_id, "n_realizations", n_realizations);

  // Set k_sum, keff, and current_batch based on whether restart file is part
  // of active cycle or inactive cycle
  restart_set_keff();
  simulation::current_batch = simulation::restart_batch;

  // Check to make sure source bank is present
  if (settings::path_sourcepoint == settings::path_statepoint &&
      !source_present) {
    fatal_error("Source bank must be contained in statepoint restart file");
  }

  // Read tallies to master. If we are using Parallel HDF5, all processes
  // need to be included in the HDF5 calls.
#ifdef PHDF5
  if (true) {
#else
  if (mpi::master) {
#endif
    // Read tally results
    load_state_point_f(file_id);
  }

  // Read source if in eigenvalue mode
  if (settings::run_mode == RUN_MODE_EIGENVALUE) {

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
  // Create type for array of 3 reals
  hsize_t dims[] {3};
  hid_t triplet = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dims);

  // Create bank datatype
  hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(struct Bank));
  H5Tinsert(banktype, "wgt", HOFFSET(Bank, wgt), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "xyz", HOFFSET(Bank, xyz), triplet);
  H5Tinsert(banktype, "uvw", HOFFSET(Bank, uvw), triplet);
  H5Tinsert(banktype, "E", HOFFSET(Bank, E), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "delayed_group", HOFFSET(Bank, delayed_group), H5T_NATIVE_INT);

  H5Tclose(triplet);
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

    std::stringstream s;
    s << settings::path_output << "source." << std::setfill('0')
      << std::setw(w) << simulation::current_batch << ".h5";
    filename_ = s.str();
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
  hsize_t count[] {static_cast<hsize_t>(simulation::work)};
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
    std::vector<Bank> temp_source {simulation::source_bank.begin(),
      simulation::source_bank.begin() + simulation::work};
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
    MPI_Send(simulation::source_bank.data(), simulation::work, mpi::bank,
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
  hsize_t dims[] {static_cast<hsize_t>(simulation::work)};
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

void write_tally_results_nr(hid_t file_id)
{
  // ==========================================================================
  // COLLECT AND WRITE GLOBAL TALLIES

  hid_t tallies_group;
  if (mpi::master) {
    // Write number of realizations
    write_dataset(file_id, "n_realizations", n_realizations);

    // Write number of global tallies
    write_dataset(file_id, "n_global_tallies", N_GLOBAL_TALLIES);

    tallies_group = open_group(file_id, "tallies");
  }

  // Get pointer to global tallies
  auto gt = global_tallies();

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

  for (int i = 1; i <= n_tallies; ++i) {
    // Skip any tallies that are not active
    bool active;
    openmc_tally_get_active(i, &active);
    if (!active) continue;

    if (mpi::master && !object_exists(file_id, "tallies_present")) {
      write_attribute(file_id, "tallies_present", 1);
    }

    // Get view of accumulated tally values
    auto results = tally_results(i);
    auto values_view = xt::view(results, xt::all(), xt::all(),
      xt::range(RESULT_SUM, RESULT_SUM_SQ + 1));

    // Make copy of tally values in contiguous array
    xt::xtensor<double, 2> values = values_view;

    if (mpi::master) {
      // Open group for tally
      int id;
      openmc_tally_get_id(i, &id);
      std::string groupname {"tally " + std::to_string(id)};
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
      xt::xtensor<double, 3> results_copy = xt::zeros_like(results);
      auto copy_view = xt::view(results_copy, xt::all(), xt::all(),
        xt::range(RESULT_SUM, RESULT_SUM_SQ + 1));
      copy_view = values;

      // Write reduced tally results to file
      auto shape = results_copy.shape();
      write_tally_results(tally_group, shape[0], shape[1], results_copy.data());

      close_group(tally_group);
    } else {
      // Receive buffer not significant at other processors
#ifdef OPENMC_MPI
      MPI_Reduce(values.data(), nullptr, values.size(), MPI_REAL8, MPI_SUM,
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
