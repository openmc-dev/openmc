#include "openmc/tallies/filter_particle.h"

#include <fmt/core.h>

#include "openmc/xml_interface.h"

namespace openmc {

void ParticleFilter::from_xml(pugi::xml_node node)
{
  auto particles = get_node_array<std::string>(node, "bins");

  // Convert to vector of ParticleType
  vector<ParticleType> types;
  for (auto& p : particles) {
    types.push_back(str_to_particle_type(p));
  }
  this->set_particles(types);
}

void ParticleFilter::set_particles(span<ParticleType> particles)
{
  // Clear existing particles
  particles_.clear();
  particles_.reserve(particles.size());

  // Set particles and number of bins
  for (auto p : particles) {
    particles_.push_back(p);
  }
  n_bins_ = particles_.size();
}

void ParticleFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  for (auto i = 0; i < particles_.size(); i++) {
    if (particles_[i] == p.type()) {
      match.bins_.push_back(i);
      match.weights_.push_back(1.0);
    }
  }
}

void ParticleFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  vector<std::string> particles;
  for (auto p : particles_) {
    particles.push_back(particle_type_to_str(p));
  }
  write_dataset(filter_group, "bins", particles);
}

std::string ParticleFilter::text_label(int bin) const
{
  const auto& p = particles_.at(bin);
  return fmt::format("Particle: {}", particle_type_to_str(p));
}

extern "C" int openmc_particle_filter_get_bins(int32_t idx, int bins[])
{
  if (int err = verify_filter(idx))
    return err;

  const auto& f = model::tally_filters[idx];
  auto pf = dynamic_cast<ParticleFilter*>(f.get());
  if (pf) {
    const auto& particles = pf->particles();
    for (int i = 0; i < particles.size(); i++) {
      bins[i] = static_cast<int>(particles[i]);
    }
  } else {
    set_errmsg("The filter at the specified index is not a ParticleFilter");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  return 0;
}

} // namespace openmc
