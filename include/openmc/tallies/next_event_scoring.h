#ifndef OPENMC_TALLIES_NEXT_EVENT_SCORING_H
#define OPENMC_TALLIES_NEXT_EVENT_SCORING_H

#include <cmath>   // for exp
#include <cstdint> // for uint64_t
#include <cstring> // for memcpy

#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/ray.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/thermal.h"

namespace openmc {

//! RNG substream offset for a detector, derived from its position.
//!
//! Keying on the position rather than on the detector's index in
//! model::active_point_detectors matters: that container is a std::set, so
//! adding a detector can renumber the ones already there. A position-derived
//! offset makes each detector's contribution depend only on where it is.
//! FNV-1a over the coordinate bits, masked to 56 bits to keep the skip-ahead
//! cheap while leaving collisions between substreams negligible.
inline int64_t point_detector_substream_offset(const Position& det)
{
  uint64_t h {14695981039346656037ULL};
  const double xyz[3] {det.x, det.y, det.z};
  for (int i = 0; i < 3; ++i) {
    uint64_t bits;
    std::memcpy(&bits, &xyz[i], sizeof(bits));
    for (int b = 0; b < 8; ++b) {
      h ^= (bits >> (8 * b)) & 0xFFULL;
      h *= 1099511628211ULL;
    }
  }
  return static_cast<int64_t>(h & 0x00FFFFFFFFFFFFFFULL);
}

//==============================================================================
//! Switches a particle onto the dedicated point-detector RNG stream for the
//! duration of next-event scoring, restoring the previous stream on exit.
//!
//! Evaluating a next-event contribution consumes random numbers (the outgoing
//! energy is sampled, and the double-valued CM->lab branch picks a root). Those
//! draws must not come from STREAM_TRACKING: otherwise defining a point tally,
//! or adding a second detector, perturbs the transport itself.
//==============================================================================

class PointDetectorStream {
public:
  explicit PointDetectorStream(Particle& p) : p_ {p}, saved_stream_ {p.stream()}
  {
    p_.stream() = STREAM_POINT_DETECTOR;
  }

  ~PointDetectorStream() { p_.stream() = saved_stream_; }

  PointDetectorStream(const PointDetectorStream&) = delete;
  PointDetectorStream& operator=(const PointDetectorStream&) = delete;

private:
  Particle& p_;
  int saved_stream_;
};

void score_point_tally_elastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, Direction v_t);

void score_point_tally_inelastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, double yield);

void score_point_tally_fission(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product);

void score_point_tally_sab(Particle& p, int i_nuclide, const ThermalData& sab,
  const NuclideMicroXS& micro);

void score_point_tally_source(SourceSite& site, int source_index);

//! Accumulate next-event contributions to every active point detector.
//!
//! \param[in] r Position of the emitting event
//! \param[in] type Type of the emitted particle
//! \param[in] time Time of the emitting event
//! \param[in] wgt Weight of the emitting particle (or source site)
//! \param[inout] seed Point-detector RNG stream, or nullptr if pdffunc
//!   consumes no random numbers
//! \param[in] pdffunc Returns the emission density per steradian toward a
//!   given direction, and sets the outgoing energy
template<typename PDF>
void score_point_tally_impl(const Position r, const ParticleType type,
  const double time, const double wgt, uint64_t* seed, PDF pdffunc)
{
  // A particle carrying no weight makes no contribution, and would divide by
  // zero in PointFilter::get_all_bins
  if (wgt == 0.0)
    return;

  // Base of this event's substream. Each detector is given its own offset from
  // it so that the contribution scored to a given detector does not depend on
  // how many other detectors are present, nor on how many random numbers the
  // preceding detectors happened to consume.
  const uint64_t seed_start = (seed != nullptr) ? *seed : 0;

  for (auto& det : model::active_point_detectors) {
    if (seed != nullptr) {
      *seed = seed_start;
      advance_prn_seed(point_detector_substream_offset(det), seed);
    }

    auto u = (det - r);
    double total_distance = u.norm();
    u /= total_distance;
    double E;
    double pdf = pdffunc(u, E);
    if (pdf == 0.0)
      continue;
    auto p = ParticleRay(r, u, type, time, E);
    p.Ray::trace(total_distance);
    // The ray left the model before reaching the detector
    if (!p.completed())
      continue;
    double mfp = p.traversal_mfp();
    double attenuation = std::exp(-mfp);

    // Carry the emitting particle's weight, and save the attenuation for point
    // filter handling
    p.wgt() = wgt;
    p.wgt_last() = wgt;
    p.wgt() *= attenuation;

    double flux = p.wgt_last() * pdf;
    score_tracklength_tally_general(p, flux, model::active_point_tallies);
  }

  // Advance the base substream by a single step so that the next event draws
  // from fresh substreams, again independently of the number of detectors.
  if (seed != nullptr) {
    *seed = seed_start;
    advance_prn_seed(1, seed);
  }
}

} // namespace openmc

#endif // OPENMC_TALLIES_NEXT_EVENT_SCORING_H
