#include "openmc/tallies/filter_musurface.h"

#include <cmath> // for abs, copysign

#include "openmc/search.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

void MuSurfaceFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Use the normal recorded for the crossing being scored. It is already a
  // unit vector expressed in the root coordinate frame, which is the frame
  // p.u() is in -- recomputing it here from the surface would give the normal
  // in the local frame of whichever universe holds the surface.
  Direction n = p.surface_normal();

  // Determine whether normal should be pointing in or out
  if (p.surface() < 0)
    n *= -1;

  // Determine cosine of angle between normal and particle direction
  double mu = p.u().dot(n);
  if (std::abs(mu) > 1.0)
    mu = std::copysign(1.0, mu);

  // Find matching bin
  if (mu >= bins_.front() && mu <= bins_.back()) {
    auto bin = lower_bound_index(bins_.begin(), bins_.end(), mu);
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
  }
}

} // namespace openmc
