#include "openmc/tallies/filter_particle.h"

#include "openmc/xml_interface.h"

namespace openmc {

void
ParticleFilter::from_xml(pugi::xml_node node)
{
  particles_ = get_node_array<int>(node, "bins");
  for (auto& p : particles_) --p;
  n_bins_ = particles_.size();
}

void
ParticleFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  for (auto i = 0; i < particles_.size(); i++) {
    if (particles_[i] == p->type) {
      //TODO: off-by-one
      match.bins_.push_back(i + 1);
      match.weights_.push_back(1.0);
    }
  }
}

void
ParticleFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", particles_);
}

std::string
ParticleFilter::text_label(int bin) const
{
  return "Particle " + std::to_string(particles_[bin]);
}

} // namespace openmc
