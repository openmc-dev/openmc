#ifndef OPENMC_RAY_H
#define OPENMC_RAY_H

#include <algorithm> // for copy

#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/particle_type.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
//! Per-flight bookkeeping, kept separate from any geometry or physics state.
//!
//! This lives in its own class so that ParticleRay can be both a Particle and
//! a ray without the two meeting at a *virtual* GeometryState base. A virtual
//! base would put a vtable indirection in front of every geometry access --
//! r(), u(), coord(), material() -- throughout the transport loop, which is
//! far too high a price to pay across all of OpenMC for one kind of ray.
//!
//! Members are public because trace_ray() drives them directly.
//==============================================================================

class RayState {
public:
  //! Accumulate a segment of flight path, splitting it between the distance
  //! travelled overall and the distance travelled inside the model.
  //! \param distance true geometric length of the segment
  void accumulate_distance(double distance);

  // Stops the ray and exits tracing when called from on_intersection
  void stop() { stop_ = true; }

  //! Whether the ray travelled the full max_distance passed to trace().
  //! False if it left the model (or hit a dead end) beforehand. This is an
  //! exact flag rather than a comparison of total_distance() against the
  //! requested distance, which is only accurate to floating-point roundoff.
  bool completed() const { return completed_; }

  //! Distance travelled inside the model, i.e. measured from the point where
  //! the ray entered the model rather than from its origin. This is the
  //! quantity reported to on_intersection().
  double traversal_distance() const { return traversal_distance_; }

  //! Distance travelled since the ray's origin, including any flight through
  //! void before the model was reached. This is what max_distance is measured
  //! against.
  double total_distance() const { return total_distance_; }

  //! Reset everything a flight accumulates.
  void reset_ray_state() { *this = RayState {}; }

  // Records how far the ray has traveled inside the model
  double traversal_distance_ {0.0};

  // Records how far the ray has traveled in total, including the leg through
  // void before it reached the model
  double total_distance_ {0.0};

  // Whether the ray has entered the model yet. Segments flown before that
  // count toward total_distance_ but not toward traversal_distance_.
  bool in_model_ {false};

  // Set when trace() consumed the whole max_distance it was given
  bool completed_ {false};

  bool stop_ {false};

  unsigned event_counter_ {0};
};

//! Trace a ray through the geometry, calling ray.on_intersection() at every
//! surface boundary. All per-flight accumulators are reset on entry, so a ray
//! may be traced more than once.
//!
//! RayT must be a GeometryState that also carries a RayState and supplies
//! on_intersection(), update_distance(double) and reset_trace_state().
//! Explicitly instantiated in ray.cpp for Ray and ParticleRay; being a
//! template rather than a virtual method is what lets the two share this
//! logic without sharing a base class.
template<typename RayT>
void trace_ray(RayT& ray, double max_distance);

// Base class that implements ray tracing logic, not necessarily through
// defined regions of the geometry but also outside of it.
class Ray : public GeometryState, public RayState {

public:
  Ray() = default;

  // Initialize from location and direction
  Ray(Position r, Direction u) { init_from_r_u(r, u); }

  //! Initialize from a known geometry state.
  explicit Ray(const GeometryState& p)
  {
    static_cast<GeometryState&>(*this) = p;
  }

  virtual ~Ray() = default;

  // Called at every surface intersection within the model
  virtual void on_intersection() = 0;

  void trace(double max_distance = INFTY);

  //! Recompute the distance to the next surface, e.g. after a direction change
  void compute_distance() { boundary() = distance_to_boundary(*this); }

  //! Accumulate a segment of flight path.
  //! \param distance true geometric length of the segment
  virtual void update_distance(double distance)
  {
    accumulate_distance(distance);
  }

  //! Reset everything trace() accumulates. Derived classes that add their own
  //! accumulators must override this and call the base version first.
  virtual void reset_trace_state()
  {
    reset_ray_state();
    boundary().reset();
  }
};

class ParticleRay : public Particle, public RayState {

public:
  ParticleRay() = default;

  //! Re-initialize an existing ray for a new flight.
  //!
  //! Everything trace() accumulates is reset, but the heap buffers that
  //! ParticleData's constructor sizes from the loaded model -- the
  //! microscopic cross section caches, the filter matches, the flux
  //! derivatives -- are left allocated so they can be reused. That is the
  //! point of resetting rather than constructing: a caller that launches many
  //! rays would otherwise allocate and zero tens of kilobytes apiece, which is
  //! not affordable anywhere near the transport loop.
  //!
  //! Note that init_from_r_u() clears material(), so the first segment of the
  //! new flight always recomputes cross sections rather than reusing whatever
  //! the previous flight happened to leave cached.
  void reset(
    Position r, Direction u, ParticleType type_, double time_, double E_)
  {
    init_from_r_u(r, u);
    reset_trace_state();
    init_physics(type_, time_, E_);
    zero_flux_derivs();
  }

  //! Construct a free-standing ray with an explicitly chosen RNG seed.
  //
  //! \param seed_id value used to stride the ray's RNG seeds. Rays that are
  //!   spawned by a particle should use the parent-particle constructor
  //!   below instead; this one is for rays with no parent.
  ParticleRay(Position r, Direction u, ParticleType type_, double time_,
    double E_, int64_t seed_id = 0)
  {
    init_from_r_u(r, u);
    init_seeds(seed_id);
    init_physics(type_, time_, E_);
  }

  //! Construct a ray emitted by an existing particle.
  //
  //! The parent's RNG seeds are *copied*, not shared. Cross section lookups
  //! along the ray (URR probability tables, temperature interpolation,
  //! S(a,b) interpolation) consume random numbers, and drawing them from the
  //! parent's live state would perturb the parent's random walk and change
  //! the results of the underlying simulation. Copying also means the ray
  //! samples the same URR probability-table realization the parent would
  //! have, which preserves the correlation of the resonance structure along
  //! the flight path.
  ParticleRay(const Particle& parent, Direction u, double E_)
  {
    // The direction differs from the parent's, so the cached coordinate
    // levels are no longer valid and the cell has to be found again. Starting
    // from r/u rather than copying the parent's GeometryState keeps this
    // correct at the cost of one exhaustive cell search.
    init_from_r_u(parent.r(), u);
    std::copy(parent.seeds(), parent.seeds() + N_STREAMS, seeds());
    stream() = STREAM_TRACKING;
    init_physics(parent.type(), parent.time(), E_);
  }

  void trace(double max_distance = INFTY);

  //! No-op: a ParticleRay has nothing to do at a surface crossing, since the
  //! optical depth is accumulated in update_distance() instead.
  void on_intersection() {}

  void update_distance(double distance);

  //! A ray that cannot be located is not a lost particle.
  //
  //! Particle::mark_as_lost() writes a restart file, increments the shared
  //! simulation::n_lost_particles counter and can abort the run outright.
  //! None of that is appropriate for a ray, which is not transported and
  //! whose failure to find a cell simply ends the trace.
  void mark_as_lost(const char* message) override;

  // Keep the std::string / std::stringstream overloads visible; declaring
  // mark_as_lost above would otherwise hide them.
  using GeometryState::mark_as_lost;

  // Records how many mean free paths the ray traveled
  double traversal_mfp() const { return traversal_mfp_; }

  //! Which entry of model::active_point_detectors this flight is aimed at, or
  //! C_NONE for a ray that is not headed for a detector. PointFilter uses it
  //! to pick the bin without having to recognise the detector by position.
  int detector_index() const { return detector_index_; }
  void set_detector_index(int i) { detector_index_ = i; }

  void reset_trace_state()
  {
    reset_ray_state();
    boundary().reset();
    traversal_mfp_ = 0.0;
    time() = time_start_;
  }

protected:
  // Records how much mean free paths the ray traveled
  double traversal_mfp_ {0.0};

  int detector_index_ {C_NONE};

private:
  //! Give the ray a deterministic, self-contained RNG state.
  void init_seeds(int64_t seed_id)
  {
    init_particle_seeds(seed_id, seeds());
    stream() = STREAM_TRACKING;
  }

  //! Set type, time and energy (and, in multigroup mode, the energy group).
  void init_physics(ParticleType type_, double time_, double E_);

  // Time the ray was launched, restored by reset_trace_state()
  double time_start_ {0.0};
};

} // namespace openmc
#endif // OPENMC_RAY_H
