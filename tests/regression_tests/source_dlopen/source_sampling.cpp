#include <iostream>
#include <memory>

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class Source : openmc::CustomSource
{
  openmc::Particle::Bank sample(uint64_t *seed)
  {
    openmc::Particle::Bank particle;
    // wgt
    particle.particle = openmc::Particle::Type::neutron;
    particle.wgt = 1.0;
    // position
    
    particle.r.x = 0.;
    particle.r.y = 0.;
    particle.r.z = 0.;
    // angle
    particle.u = {1.0, 0.0, 0.0};
    particle.E = 14.08e6;
    particle.delayed_group = 0;
    return particle;    
  }
};

// A function to create a unique pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<Source> openmc_create_source(std::string parameters)
{
  return std::make_unique<Source>();
}
