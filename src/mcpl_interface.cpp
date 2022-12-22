#include "openmc/mcpl_interface.h"

#include "openmc/error.h"

#ifdef OPENMC_MCPL
#include <mcpl.h>
#endif

namespace openmc {

#ifdef OPENMC_MCPL
const bool MCPL_ENABLED = true;
#else
const bool MCPL_ENABLED = false;
#endif

#ifdef OPENMC_MCPL
SourceSite mcpl_particle_to_site(const mcpl_particle* particle)
{
  SourceSite site;

  switch (particle->pdgcode) {
  case 2112:
    site.particle = ParticleType::neutron;
    break;
  case 22:
    site.particle = ParticleType::photon;
    break;
  case 11:
    site.particle = ParticleType::electron;
    break;
  case -11:
    site.particle = ParticleType::positron;
    break;
  }

  // Copy position and direction
  site.r.x = mcpl_particle->position[0];
  site.r.y = mcpl_particle->position[1];
  site.r.z = mcpl_particle->position[2];
  site.u.x = mcpl_particle->direction[0];
  site.u.y = mcpl_particle->direction[1];
  site.u.z = mcpl_particle->direction[2];

  // mcpl stores kinetic energy in MeV
  site.E = mcpl_particle->ekin * 1e6;
  // mcpl stores time in ms
  site.time = mcpl_particle->time * 1e-3;
  site.wgt = mcpl_particle->weight;

  return site;
}
#endif

vector<SourceSite> mcpl_source_sites(std::string path)
{
  vector<SourceSite> sites;

#ifdef OPENMC_MCPL
  size_t n_sites = mcpl_hdr_nparticles(mcpl_file);

  for (int i = 0; i < n_sites; i++) {
    SourceSite site;

    // Extract particle from mcpl-file, checking if it is a neutron, photon,
    // electron, or positron. Otherwise skip.
    const mcpl_particle_t* particle;
    int pdg = 0;
    while (pdg != 2112 && pdg != 22 && pdg != 11 && pdg != -11) {
      particle = mcpl_read(mcpl_file);
      pdg = mcpl_particle->pdgcode;
    }

    // Convert to source site and add to vector
    sites.push_back(mcpl_particle_to_site(particle));
  }

  // Check that some sites were read
  if (sites_.empty()) {
    fatal_error("MCPL file contained no neutron, photon, electron, or positron "
                "source particles.");
  }

  mcpl_close_file(mcpl_file);
#else
  fatal_error(
    "Your build of OpenMC does not support reading MCPL source files.");
#endif

  return sites;
}

} // namespace openmc
