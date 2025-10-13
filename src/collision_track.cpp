#include "openmc/collision_track.h"

#include <algorithm>

#include <fmt/format.h>

#include "openmc/bank.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mcpl_interface.h"
#include "openmc/message_passing.h"
#include "openmc/output.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc::collision_track {

namespace {

hid_t h5_collision_track_banktype()
{
  hid_t postype = H5Tcreate(H5T_COMPOUND, sizeof(Position));
  H5Tinsert(postype, "x", HOFFSET(Position, x), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "y", HOFFSET(Position, y), H5T_NATIVE_DOUBLE);
  H5Tinsert(postype, "z", HOFFSET(Position, z), H5T_NATIVE_DOUBLE);

  hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(CollisionTrackSite));

  H5Tinsert(banktype, "r", HOFFSET(CollisionTrackSite, r), postype);
  H5Tinsert(banktype, "u", HOFFSET(CollisionTrackSite, u), postype);
  H5Tinsert(banktype, "E", HOFFSET(CollisionTrackSite, E), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "dE", HOFFSET(CollisionTrackSite, dE), H5T_NATIVE_DOUBLE);
  H5Tinsert(
    banktype, "time", HOFFSET(CollisionTrackSite, time), H5T_NATIVE_DOUBLE);
  H5Tinsert(
    banktype, "wgt", HOFFSET(CollisionTrackSite, wgt), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "event_mt", HOFFSET(CollisionTrackSite, event_mt),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "delayed_group",
    HOFFSET(CollisionTrackSite, delayed_group), H5T_NATIVE_INT);
  H5Tinsert(
    banktype, "cell_id", HOFFSET(CollisionTrackSite, cell_id), H5T_NATIVE_INT);
  H5Tinsert(banktype, "nuclide_id", HOFFSET(CollisionTrackSite, nuclide_id),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "material_id", HOFFSET(CollisionTrackSite, material_id),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "universe_id", HOFFSET(CollisionTrackSite, universe_id),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "n_collision", HOFFSET(CollisionTrackSite, n_collision),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "particle", HOFFSET(CollisionTrackSite, particle),
    H5T_NATIVE_INT);
  H5Tinsert(banktype, "parent_id", HOFFSET(CollisionTrackSite, parent_id),
    H5T_NATIVE_LONG);
  H5Tinsert(banktype, "progeny_id", HOFFSET(CollisionTrackSite, progeny_id),
    H5T_NATIVE_LONG);
  H5Tclose(postype);
  return banktype;
}

void write_collision_track_bank(hid_t group_id,
  openmc::span<CollisionTrackSite> collision_track_bank,
  const openmc::vector<int64_t>& bank_index)
{
  hid_t banktype = h5_collision_track_banktype();

  // Set total and individual process dataspace sizes for source bank
  int64_t dims_size = bank_index.back();
  int64_t count_size = bank_index[mpi::rank + 1] - bank_index[mpi::rank];

#ifdef PHDF5
  // Set size of total dataspace for all procs and rank
  hsize_t dims[] {static_cast<hsize_t>(dims_size)};
  hid_t dspace = H5Screate_simple(1, dims, nullptr);
  hid_t dset = H5Dcreate(group_id, "collision_track_bank", banktype, dspace,
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create another data space but for each proc individually
  hsize_t count[] {static_cast<hsize_t>(count_size)};
  hid_t memspace = H5Screate_simple(1, count, nullptr);

  // Select hyperslab for this dataspace
  hsize_t start[] {static_cast<hsize_t>(bank_index[mpi::rank])};
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Set up the property list for parallel writing
  hid_t plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

  // Write data to file in parallel
  H5Dwrite(
    dset, banktype, memspace, dspace, plist, collision_track_bank.data());

  // Free resources
  H5Sclose(dspace);
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Pclose(plist);

#else

  if (mpi::master) {
    // Create dataset big enough to hold all collisions
    hsize_t dims[] {static_cast<hsize_t>(dims_size)};
    hid_t dspace = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate(group_id, "collision_track_bank", banktype, dspace,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Save collision track events since the array is overwritten below
#ifdef OPENMC_MPI
    openmc::vector<CollisionTrackSite> temp_collision_track {
      collision_track_bank.begin(), collision_track_bank.end()};
#endif

    for (int i = 0; i < mpi::n_procs; ++i) {
      // Create memory space
      hsize_t count[] {static_cast<hsize_t>(bank_index[i + 1] - bank_index[i])};
      hid_t memspace = H5Screate_simple(1, count, nullptr);

#ifdef OPENMC_MPI
      // Receive collision track sites from other processes
      if (i > 0)
        MPI_Recv(collision_track_bank.data(), count[0],
          mpi::collision_track_site, i, i, mpi::intracomm, MPI_STATUS_IGNORE);
#endif

      // Select hyperslab for this dataspace
      dspace = H5Dget_space(dset);
      hsize_t start[] {static_cast<hsize_t>(bank_index[i])};
      H5Sselect_hyperslab(
        dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

      // Write data to hyperslab
      H5Dwrite(dset, banktype, memspace, dspace, H5P_DEFAULT,
        collision_track_bank.data());

      H5Sclose(memspace);
      H5Sclose(dspace);
    }

    // Close all ids
    H5Dclose(dset);

#ifdef OPENMC_MPI
    // Restore state of collision track bank
    std::copy(temp_collision_track.begin(), temp_collision_track.end(),
      collision_track_bank.begin());
#endif
  } else {
#ifdef OPENMC_MPI
    MPI_Send(collision_track_bank.data(), count_size, mpi::collision_track_site,
      0, mpi::rank, mpi::intracomm);
#endif
  }
#endif

  H5Tclose(banktype);
}

void write_h5_collision_track(const char* filename,
  openmc::span<CollisionTrackSite> collision_track_bank,
  const openmc::vector<int64_t>& bank_index)
{
#ifdef PHDF5
  bool parallel = true;
#else
  bool parallel = false;
#endif

  if (!filename)
    fatal_error("write_h5_collision_track filename needs a nonempty name.");

  std::string filename_(filename);
  const auto extension = get_file_extension(filename_);
  if (extension.empty()) {
    filename_.append(".h5");
  } else if (extension != "h5") {
    warning("write_h5_collision_track was passed a file extension differing "
            "from .h5, but an hdf5 file will be written.");
  }

  hid_t file_id;
  if (mpi::master || parallel) {
    file_id = file_open(filename_.c_str(), 'w', true);
    write_attribute(file_id, "filetype", "source");
  }

  write_collision_track_bank(file_id, collision_track_bank, bank_index);

  if (mpi::master || parallel)
    file_close(file_id);
}

} // namespace

RuntimeState state {};

void reset_runtime()
{
  state.current_file = 1;
}

bool should_record_event(int id_cell, int mt_event, const std::string& nuclide,
  int id_universe, int id_material, double energy_loss)
{
  auto matches_filter = [](const auto& filter_set, const auto& value) {
    return filter_set.empty() || filter_set.count(value) > 0;
  };

  const auto& cfg = settings::collision_track_config;
  return simulation::current_batch > settings::n_inactive &&
         !simulation::collision_track_bank.full() &&
         matches_filter(cfg.cell_ids, id_cell) &&
         matches_filter(cfg.mt_numbers, mt_event) &&
         matches_filter(cfg.universe_ids, id_universe) &&
         matches_filter(cfg.material_ids, id_material) &&
         matches_filter(cfg.nuclides, nuclide) &&
         (cfg.deposited_energy_threshold == 0 ||
           cfg.deposited_energy_threshold < energy_loss);
}

void reserve_bank_capacity()
{
  simulation::collision_track_bank.reserve(
    settings::collision_track_config.max_collisions);
}

void flush_bank(bool last_batch)
{
  if (!settings::collision_track)
    return;

  const auto& cfg = settings::collision_track_config;
  if (state.current_file > cfg.max_files)
    return;

  if (!simulation::collision_track_bank.full() && !last_batch)
    return;

  auto size = simulation::collision_track_bank.size();
  if (size == 0 && !last_batch)
    return;

  auto collision_track_work_index = mpi::calculate_parallel_index_vector(size);
  openmc::span<CollisionTrackSite> collisiontrackbankspan(
    simulation::collision_track_bank.begin(), size);

  std::string ext = cfg.mcpl_write ? "mcpl" : "h5";
  auto filename = fmt::format(
    "{}collision_track.{}.{}", settings::path_output, state.current_file, ext);

  if (cfg.max_files == 1 || (state.current_file == 1 && last_batch)) {
    filename = settings::path_output + "collision_track." + ext;
    write_message(filename + " file with {} recorded collisions ...", size, 2);
  } else {
    write_message(
      "Creating collision_track.{}.{} file with {} recorded collisions ...",
      state.current_file, ext, size, 4);
  }

  if (cfg.mcpl_write) {
    write_mcpl_collision_track(
      filename.c_str(), collisiontrackbankspan, collision_track_work_index);
  } else {
    write_h5_collision_track(
      filename.c_str(), collisiontrackbankspan, collision_track_work_index);
  }

  simulation::collision_track_bank.clear();
  if (!last_batch && cfg.max_files >= 1) {
    reserve_bank_capacity();
  }
  ++state.current_file;
}

void reset_config()
{
  settings::collision_track_config = settings::CollisionTrackConfig {};
  reset_runtime();
}

} // namespace openmc::collision_track
