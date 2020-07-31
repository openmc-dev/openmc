#include <cmath> // for M_PI

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"
#include "pugixml.hpp"

class SerialisedSource {
  protected:
    double radius_;
    double energy_;

    // Protect the constructor so that the class can only be created by serialisation.
    SerialisedSource(double radius, double energy) {
      radius_ = radius;
      energy_ = energy;
    }

  public:
    // Getters for the values that we want to use in sampling.
    double radius() { return radius_; }
    double energy() { return energy_; }

    // The deserialisation routine populates the constructor from well defined elements
    // in the input XML document.
    // Note that the source will have already been read from file, so what will be passed
    // in here is a string-like serialized value (not the path to the serialized value).
    static SerialisedSource from_xml(char* serialised_source) {
      pugi::xml_document doc;
      doc.load_string(serialised_source);
      pugi::xml_node root_node = doc.root().child("Source");
      double radius = root_node.child("Radius").text().as_double();
      double energy = root_node.child("Energy").text().as_double();
      return SerialisedSource(radius, energy);
    }
};

// you must have external C linkage here otherwise
// dlopen will not find the file
extern "C" openmc::Particle::Bank sample_source(uint64_t* seed, char* serialised_source)
{
  SerialisedSource source = SerialisedSource::from_xml(serialised_source);

  openmc::Particle::Bank particle;
  // wgt
  particle.particle = openmc::Particle::Type::neutron;
  particle.wgt = 1.0;
  // position
  double angle = 2. * M_PI * openmc::prn(seed);
  double radius = source.radius();
  particle.r.x = radius * std::cos(angle);
  particle.r.y = radius * std::sin(angle);
  particle.r.z = 0.0;
  // angle
  particle.u = {1.0, 0.0, 0.0};
  particle.E = source.energy();
  particle.delayed_group = 0;

  return particle;
}
