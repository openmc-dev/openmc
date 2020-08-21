#include "openmc/source.h"
#include "openmc/particle.h"

class Source : public openmc::CustomSource
{
  public:
    double energy;

    Source(double energy)
    {
      this->energy = energy;
    }

    // Samples from an instance of this class.
    openmc::Particle::Bank sample_source(uint64_t* seed)
    {
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

// A function to create a unique pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<Source> openmc_create_source(const char* parameter)
{
  return std::unique_ptr<Source> (new Source(atof(parameter)));
}
