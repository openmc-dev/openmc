#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"
#include "pugixml.hpp"

class SerializedSource {
  protected:
    double energy_;

    // Protect the constructor so that the class can only be created by serialisation.
    SerializedSource(double energy) {
      energy_ = energy;
    }

  public:
    // Getters for the values that we want to use in sampling.
    double energy() { return energy_; }

    static SerializedSource from_xml(char* serialized_source) {
      pugi::xml_document doc;
      doc.load_string(serialized_source);
      pugi::xml_node root_node = doc.root().child("Source");
      double energy = root_node.child("Energy").text().as_double();
      return SerializedSource(energy);
    }
};

// you must have external C linkage here otherwise
// dlopen will not find the file
extern "C" openmc::Particle::Bank sample_source(uint64_t* seed, char* serialized_source) {
  SerializedSource source = SerializedSource::from_xml(serialized_source);

  openmc::Particle::Bank particle;
  // wgt
  particle.particle = openmc::Particle::Type::neutron;
  particle.wgt = 1.0;
  // position
  particle.r.x = 0.0;
  particle.r.y = 0.0;
  particle.r.z = 0.0;
  // angle
  particle.u = {1.0, 0.0, 0.0};
  particle.E = source.energy();
  particle.delayed_group = 0;

  return particle;
}
