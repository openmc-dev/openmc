#ifndef OPENMC_RAY_H
#define OPENMC_RAY_H

#include "openmc/particle_data.h"
#include "openmc/position.h"

namespace openmc {

// Base class that implements ray tracing logic, not necessarily through
// defined regions of the geometry but also outside of it.
class Ray : public GeometryState {

public:
  // Initialize from location and direction
  Ray(Position r, Direction u) { init_from_r_u(r, u); }

  // Initialize from known geometry state
  Ray(const GeometryState& p) : GeometryState(p) {}

  // Called at every surface intersection within the model
  virtual void on_intersection() = 0;

  /*
   * Traces the ray through the geometry, calling on_intersection
   * at every surface boundary.
   */
  void trace();

  // Stops the ray and exits tracing when called from on_intersection
  void stop() { stop_ = true; }

  // Sets the dist_ variable
  void compute_distance();

protected:
  // Records how far the ray has traveled
  double traversal_distance_ {0.0};

private:
  // Max intersections before we assume ray tracing is caught in an infinite
  // loop:
  static const int MAX_INTERSECTIONS = 1000000;

  bool stop_ {false};

  unsigned event_counter_ {0};
};

} // namespace openmc
#endif // OPENMC_RAY_H
