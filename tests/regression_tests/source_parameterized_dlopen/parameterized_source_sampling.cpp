#include "openmc/source.h"
#include "openmc/particle.h"

class ParameterizedSource : public openmc::CustomSource {
  public:
    double energy;

    ParameterizedSource(double energy) {
      this->energy = energy;
    }

    // Samples from an instance of this class.
    openmc::Particle::Bank sample_source(uint64_t* seed) {
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
      particle.E = this->energy;
      particle.delayed_group = 0;

      return particle;
    }
};

// you must have external C linkage here otherwise
// dlopen will not find the file
extern "C" ParameterizedSource* openmc_create_source(const char* parameter) {
  return new ParameterizedSource(atof(parameter));
}
