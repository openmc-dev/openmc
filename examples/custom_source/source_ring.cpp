#include <cmath> // for M_PI
#include <memory> // for unique_ptr

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class RingSource : public openmc::Source
{
  openmc::SourceSite sample(uint64_t* seed) const
  {
    openmc::SourceSite particle;
    // particle type
    particle.particle = openmc::ParticleType::neutron;
    // position
    double angle = 2.0 * M_PI * openmc::prn(seed);
    double radius = 3.0;
    particle.r.x = radius * std::cos(angle);
    particle.r.y = radius * std::sin(angle);
    particle.r.z = 0.0;
    // angle
    particle.u = {1.0, 0.0, 0.0};
    particle.E = 14.08e6;
    return particle;
  }
};

// A function to create a unique pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<RingSource> openmc_create_source(std::string parameters)
{
  return std::make_unique<RingSource>();
}
