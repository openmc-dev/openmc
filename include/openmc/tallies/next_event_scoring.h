#ifndef OPENMC_TALLIES_NEXT_EVENT_SCORING_H
#define OPENMC_TALLIES_NEXT_EVENT_SCORING_H

#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/ray.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/thermal.h"

namespace openmc {

void score_point_tally_elastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, Direction v_t);

void score_point_tally_inelastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, double yield);

void score_point_tally_fission(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product);

void score_point_tally_sab(Particle& p, int i_nuclide, const ThermalData& sab,
  const NuclideMicroXS& micro);

void score_point_tally_source(SourceSite& site, int source_index);

template<typename PDF>
void score_point_tally_impl(
  const Position r, const ParticleType type, const double time, PDF pdffunc)
{
  for (auto& det : model::active_point_detectors) {
    auto u = (det - r);
    double total_distance = u.norm();
    u /= total_distance;
    double E;
    double pdf = pdffunc(u, E);
    if (pdf == 0.0)
      continue;
    auto p = ParticleRay(r, u, type, time, E);
    p.Ray::trace(total_distance);
    double distance = p.traversal_distance();
    if (distance < total_distance)
      continue;
    double mfp = p.traversal_mfp();
    double attenuation = std::exp(-mfp);

    // Save the attenuation for point filter handling
    p.wgt_last() = p.wgt();
    p.wgt() *= attenuation;

    double flux = p.wgt_last() * pdf;
    score_tracklength_tally_general(p, flux, model::active_point_tallies);
  }
}

} // namespace openmc

#endif // OPENMC_TALLIES_NEXT_EVENT_SCORING_H
