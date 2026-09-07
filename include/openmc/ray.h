#ifndef OPENMC_RAY_H
#define OPENMC_RAY_H

#include <algorithm> // for copy

#include "openmc/particle.h"
#include "openmc/particle_type.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"

namespace openmc {

// Base class that implements ray tracing logic, not necessarily through
// defined regions of the geometry but also outside of it.
class Ray : virtual public GeometryState {

public:
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

  /*
   * Traces the ray through the geometry, calling on_intersection
   * at every surface boundary. All per-trace accumulators are reset on
   * entry, so a Ray may be traced more than once.
   */
  void trace(double max_distance = INFTY);

  // Stops the ray and exits tracing when called from on_intersection
  void stop() { stop_ = true; }

  // Sets the dist_ variable
  void compute_distance();

  //! Accumulate a segment of flight path.
  //! \param distance true geometric length of the segment
  virtual void update_distance(double distance);

  //! Whether the ray travelled the full max_distance passed to trace().
  //! False if it left the model (or hit a dead end) beforehand. This is an
  //! exact flag rather than a comparison of traversal_distance() against the
  //! requested distance, which is only accurate to floating-point roundoff.
  bool completed() const { return completed_; }

  // Records how far the ray has traveled, including any flight through void
  // before it reached the model
  double traversal_distance() const { return traversal_distance_; }

protected:
  //! Reset everything trace() accumulates. Derived classes that add their own
  //! accumulators must override this and call the base version first.
  virtual void reset_trace_state()
  {
    traversal_distance_ = 0.0;
    completed_ = false;
    stop_ = false;
    event_counter_ = 0;
    boundary().reset();
  }

  // Records how far the ray has traveled
  double traversal_distance_ {0.0};

  // Set when trace() consumed the whole max_distance it was given
  bool completed_ {false};

private:
  // Max intersections before we assume ray tracing is caught in an infinite
  // loop:
  static constexpr int MAX_INTERSECTIONS = 1000000;

  bool stop_ {false};

  unsigned event_counter_ {0};
};

class ParticleRay : public Ray, public Particle {

public:
  //! Construct a free-standing ray with an explicitly chosen RNG seed.
  //
  //! \param seed_id value used to stride the ray's RNG seeds. Rays that are
  //!   spawned by a particle should use the parent-particle constructor
  //!   below instead; this one is for rays with no parent.
  ParticleRay(Position r, Direction u, ParticleType type_, double time_,
    double E_, int64_t seed_id = 0)
    : Ray(r, u)
  {
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
    // The direction differs from the parent's, so the cached coordinate
    // levels are no longer valid and the cell has to be found again. Starting
    // from r/u rather than copying the parent's GeometryState keeps this
    // correct at the cost of one exhaustive cell search.
    : Ray(parent.r(), u)
  {
    std::copy(parent.seeds(), parent.seeds() + N_STREAMS, seeds());
    stream() = STREAM_TRACKING;
    init_physics(parent.type(), parent.time(), E_);
  }

  void on_intersection() override;

  void update_distance(double distance) override;

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

protected:
  void reset_trace_state() override
  {
    Ray::reset_trace_state();
    traversal_mfp_ = 0.0;
    time() = time_start_;
  }

  // Records how much mean free paths the ray traveled
  double traversal_mfp_ {0.0};

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
