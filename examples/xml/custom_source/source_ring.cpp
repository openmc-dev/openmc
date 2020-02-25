#include <cmath> // for M_PI

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

// you must have external C linkage here otherwise
// dlopen will not find the file
extern "C" openmc::Particle::Bank sample_source(uint64_t* seed)
{
  openmc::Particle::Bank particle;
  // wgt
  particle.particle = openmc::Particle::Type::neutron;
  particle.wgt = 1.0;
  // position
  double angle = 2. * M_PI * openmc::prn(seed);
  double radius = 3.0;
  particle.r.x = radius * std::cos(angle);
  particle.r.y = radius * std::sin(angle);
  particle.r.z = 0.0;
  // angle
  particle.u = {1.0, 0.0, 0.0};
  particle.E = 14.08e6;
  particle.delayed_group = 0;
  return particle;
}
