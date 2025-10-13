#ifndef OPENMC_COLLISION_TRACK_H
#define OPENMC_COLLISION_TRACK_H

#include "openmc/particle_data.h"
#include "openmc/span.h"
#include "openmc/vector.h"

#include <cstdint>
#include <string>

namespace openmc::collision_track {

struct RuntimeState {
  int current_file {1};
};

extern RuntimeState state;

void reset_runtime();

bool should_record_event(int id_cell, int mt_event, const std::string& nuclide,
  int id_universe, int id_material, double energy_loss);

void reserve_bank_capacity();

void flush_bank(bool last_batch);

void reset_config();

} // namespace openmc::collision_track

#endif // OPENMC_COLLISION_TRACK_H
