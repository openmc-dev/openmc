#include <iostream>
#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

// you must have external C linkage here otherwise 
// dlopen will not find the file
extern "C" openmc::Particle::Bank sample_source(uint64_t *seed) {
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
