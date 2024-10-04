#include "openmc/tallies/filter_musurface.h"

#include <cmath> // for abs, copysign

#include "openmc/search.h"
#include "openmc/surface.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

void MuSurfaceFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Get surface normal (and make sure it is a unit vector)
  const auto surf {model::surfaces[std::abs(p.surface()) - 1].get()};
  auto n = surf->normal(p.r());
  n /= n.norm();

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
