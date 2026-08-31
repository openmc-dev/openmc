//! \file pulse_height.h
//! \brief Deferred, per-history aggregation of pulse-height results
//!
//! A pulse-height tally scores one count per source history, in the bin
//! containing the total energy that history's entire particle tree deposited in
//! a given cell. In the default transport modes a whole tree is carried by a
//! single Particle object, so Particle::pht_storage() already holds the
//! per-history total by the time event_death() runs and can be scored directly.
//!
//! Under the shared secondary bank each secondary generation is transported as
//! a fresh set of Particle objects, redistributed across MPI ranks between
//! generations. A history's deposition is therefore spread over many Particle
//! objects on potentially many ranks. This module collects those fragments,
//! keyed by the root index carried on every SourceSite, and scores them once
//! per history after the generation loop has drained.

#ifndef OPENMC_TALLIES_PULSE_HEIGHT_H
#define OPENMC_TALLIES_PULSE_HEIGHT_H

#include <cstdint>

#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! One flushed pulse-height fragment, tagged with the history it belongs to.
//==============================================================================

struct PulseHeightContribution {
  int64_t root_index;    //!< index of the primary at the root of the tree
  vector<double> energy; //!< per-cell energy, indexed as pulse_height_cells
};

namespace simulation {

//! Per-thread staging buffers, merged in finalize_pulse_height_tallies().
extern vector<vector<PulseHeightContribution>> pht_thread_buffers;

} // namespace simulation

//! Allocate the per-thread staging buffers. Called from initialize_simulation()
//! when pulse-height tallies and the shared secondary bank are both active.
void init_pulse_height_buffers();

//! Release the staging buffers and the phase-1 partition snapshot.
void free_memory_pulse_height();

//! Stage one Particle's contribution to its history's pulse height.
//
//! Thread-safe by construction: each thread appends only to its own buffer.
//! Contributions that are identically zero in every cell are dropped; histories
//! that deposit nothing are recovered in finalize_pulse_height_tallies() by
//! iterating over the full root range rather than over staged entries.
//
//! \param root_index index of the primary at the root of this particle's tree
//! \param pht per-cell energy deposited by this particle alone
void stage_pulse_height(int64_t root_index, const vector<double>& pht);

//! Aggregate staged contributions by history and score them.
//
//! Sends each contribution to the rank that owns its root according to
//! simulation::phase1_work_index, sums per (history, cell), and scores every
//! owned history including those with no deposition. Must be called after the
//! last secondary generation has been transported and before tally results are
//! accumulated for the batch.
void finalize_pulse_height_tallies();

} // namespace openmc

#endif // OPENMC_TALLIES_PULSE_HEIGHT_H
