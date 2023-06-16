#include <iostream>
#include <memory>

#include "openmc/particle_data.h"
#include "openmc/random_lcg.h"
#include "openmc/source.h"

class CustomSource : public openmc::Source {
  openmc::SourceSite sample(uint64_t* seed) const
  {
    openmc::SourceSite particle;
    // wgt
    particle.particle = openmc::ParticleType::neutron;
    particle.wgt = 1.0;
    // position

    particle.r.x = 0.;
    particle.r.y = 0.;
    particle.r.z = 0.;
    // angle
    particle.u = {1.0, 0.0, 0.0};
    particle.E = 1.00e3;
    particle.delayed_group = 0;
    return particle;
  }
};

// A function to create a unique pointer to an instance of this class when
// generated via a plugin call using dlopen/dlsym. You must have external C
// linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<CustomSource> openmc_create_source(
  std::string parameters)
{
  return std::make_unique<CustomSource>();
}
