#include "openmc/collision_track.h"

#include <algorithm>
#include <string>

#include <fmt/format.h>

#include "openmc/bank.h"
#include "openmc/bank_io.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mcpl_interface.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/universe.h"

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {

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
    H5T_NATIVE_INT64);
  H5Tinsert(banktype, "progeny_id", HOFFSET(CollisionTrackSite, progeny_id),
    H5T_NATIVE_INT64);
  H5Tclose(postype);
  return banktype;
}

void write_collision_track_bank(hid_t group_id,
  openmc::span<CollisionTrackSite> collision_track_bank,
  const openmc::vector<int64_t>& bank_index)
{
  hid_t banktype = h5_collision_track_banktype();
#ifdef OPENMC_MPI
  write_bank_dataset("collision_track_bank", group_id, collision_track_bank,
    bank_index, banktype, mpi::collision_track_site);
#else
  write_bank_dataset("collision_track_bank", group_id, collision_track_bank,
    bank_index, banktype);
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

    // Write filetype and version info
    write_attribute(file_id, "filetype", "collision_track");
    write_attribute(file_id, "version", VERSION_COLLISION_TRACK);
  }

  write_collision_track_bank(file_id, collision_track_bank, bank_index);

  if (mpi::master || parallel)
    file_close(file_id);
}

} // namespace

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

void collision_track_reserve_bank()
{
  simulation::collision_track_bank.reserve(
    settings::collision_track_config.max_collisions);
}

void collision_track_flush_bank()
{
  const auto& cfg = settings::collision_track_config;
  if (simulation::ct_current_file > cfg.max_files)
    return;

  bool last_batch = (simulation::current_batch == settings::n_batches);
  if (!simulation::collision_track_bank.full() && !last_batch)
    return;

  auto size = simulation::collision_track_bank.size();
  if (size == 0 && !last_batch)
    return;

  auto collision_track_work_index = mpi::calculate_parallel_index_vector(size);
  openmc::span<CollisionTrackSite> collisiontrackbankspan(
    simulation::collision_track_bank.begin(), size);

  std::string ext = cfg.mcpl_write ? "mcpl" : "h5";
  auto filename = fmt::format("{}collision_track.{}.{}", settings::path_output,
    simulation::ct_current_file, ext);

  if (cfg.max_files == 1 || (simulation::ct_current_file == 1 && last_batch)) {
    filename = settings::path_output + "collision_track." + ext;
  }
  write_message("Creating {}...", filename, 4);

  if (cfg.mcpl_write) {
    write_mcpl_collision_track(
      filename.c_str(), collisiontrackbankspan, collision_track_work_index);
  } else {
    write_h5_collision_track(
      filename.c_str(), collisiontrackbankspan, collision_track_work_index);
  }

  simulation::collision_track_bank.clear();
  if (!last_batch && cfg.max_files >= 1) {
    collision_track_reserve_bank();
  }
  ++simulation::ct_current_file;
}

void collision_track_record(Particle& particle)
{
  int cell_index = particle.lowest_coord().cell();
  if (cell_index == C_NONE)
    return;

  int cell_id = model::cells[cell_index]->id_;
  const auto* nuclide_ptr = data::nuclides[particle.event_nuclide()].get();
  std::string nuclide = nuclide_ptr->name_;
  int universe_id = model::universes[particle.lowest_coord().universe()]->id_;
  double delta_E = particle.E_last() - particle.E();
  int material_index = particle.material();
  if (material_index == C_NONE)
    return;

  int material_id = model::materials[material_index]->id_;

  if (!should_record_event(cell_id, particle.event_mt(), nuclide, universe_id,
        material_id, delta_E))
    return;

  CollisionTrackSite site;
  site.r = particle.r();
  site.u = particle.u();
  site.E = particle.E_last();
  site.dE = delta_E;
  site.time = particle.time();
  site.wgt = particle.wgt();
  site.event_mt = particle.event_mt();
  site.delayed_group = particle.delayed_group();
  site.cell_id = cell_id;
  site.nuclide_id =
    10000 * nuclide_ptr->Z_ + 10 * nuclide_ptr->A_ + nuclide_ptr->metastable_;
  site.material_id = material_id;
  site.universe_id = universe_id;
  site.n_collision = particle.n_collision();
  site.particle = particle.type();
  site.parent_id = particle.id();
  site.progeny_id = particle.n_progeny();
  simulation::collision_track_bank.thread_safe_append(site);
}

} // namespace openmc
