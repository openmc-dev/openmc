#include "openmc/tallies/filter_particle.h"

#include "openmc/xml_interface.h"

namespace openmc {

void
ParticleFilter::from_xml(pugi::xml_node node)
{
  auto particles = get_node_array<std::string>(node, "bins");

  // Convert to vector of Particle::Type
  std::vector<Particle::Type> types;
  for (auto& p : particles) {
    if (p == "neutron") {
      types.push_back(Particle::Type::neutron);
    } else if (p == "photon") {
      types.push_back(Particle::Type::photon);
    } else if (p == "electron") {
      types.push_back(Particle::Type::electron);
    } else if (p == "positron") {
      types.push_back(Particle::Type::positron);
    }
  }
  this->set_particles(types);
}

void
ParticleFilter::set_particles(gsl::span<Particle::Type> particles)
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

void
ParticleFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  for (auto i = 0; i < particles_.size(); i++) {
    if (particles_[i] == p->type_) {
      match.bins_.push_back(i);
      match.weights_.push_back(1.0);
    }
  }
}

void
ParticleFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  std::vector<std::string> particles;
  for (auto p : particles_) {
    switch (p) {
    case Particle::Type::neutron:
      particles.push_back("neutron");
      break;
    case Particle::Type::photon:
      particles.push_back("photon");
      break;
    case Particle::Type::electron:
      particles.push_back("electron");
      break;
    case Particle::Type::positron:
      particles.push_back("positron");
      break;
    }
  }
  write_dataset(filter_group, "bins", particles);
}

std::string
ParticleFilter::text_label(int bin) const
{
  switch (particles_[bin]) {
  case Particle::Type::neutron:
    return "Particle: neutron";
  case Particle::Type::photon:
    return "Particle: photon";
  case Particle::Type::electron:
    return "Particle: electron";
  case Particle::Type::positron:
    return "Particle: positron";
  }
  UNREACHABLE();
}

} // namespace openmc
