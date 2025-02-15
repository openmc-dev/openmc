#include "openmc/mcpl_interface.h"

#include "openmc/bank.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/message_passing.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/vector.h"

#include <fmt/core.h>

#ifdef OPENMC_MCPL
#include <mcpl.h>
#endif

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

#ifdef OPENMC_MCPL
const bool MCPL_ENABLED = true;
#else
const bool MCPL_ENABLED = false;
#endif

//==============================================================================
// Functions
//==============================================================================

#ifdef OPENMC_MCPL
SourceSite mcpl_particle_to_site(const mcpl_particle_t* particle)
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
  site.r.x = particle->position[0];
  site.r.y = particle->position[1];
  site.r.z = particle->position[2];
  site.u.x = particle->direction[0];
  site.u.y = particle->direction[1];
  site.u.z = particle->direction[2];

  // MCPL stores kinetic energy in [MeV], time in [ms]
  site.E = particle->ekin * 1e6;
  site.time = particle->time * 1e-3;
  site.wgt = particle->weight;

  return site;
}
#endif

//==============================================================================

vector<SourceSite> mcpl_source_sites(std::string path)
{
  vector<SourceSite> sites;

#ifdef OPENMC_MCPL
  // Open MCPL file and determine number of particles
  auto mcpl_file = mcpl_open_file(path.c_str());
  size_t n_sites = mcpl_hdr_nparticles(mcpl_file);

  for (int i = 0; i < n_sites; i++) {
    // Extract particle from mcpl-file, checking if it is a neutron, photon,
    // electron, or positron. Otherwise skip.
    const mcpl_particle_t* particle;
    int pdg = 0;
    while (pdg != 2112 && pdg != 22 && pdg != 11 && pdg != -11) {
      particle = mcpl_read(mcpl_file);
      pdg = particle->pdgcode;
    }

    // Convert to source site and add to vector
    sites.push_back(mcpl_particle_to_site(particle));
  }

  // Check that some sites were read
  if (sites.empty()) {
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

//==============================================================================

#ifdef OPENMC_MCPL
void write_mcpl_source_bank(mcpl_outfile_t file_id,
  span<SourceSite> source_bank, const vector<int64_t>& bank_index)
{
  int64_t dims_size = settings::n_particles;
  int64_t count_size = simulation::work_per_rank;

  if (mpi::master) {
    // Particles are writeen to disk from the master node only

    // Save source bank sites since the array is overwritten below
#ifdef OPENMC_MPI
    vector<SourceSite> temp_source {source_bank.begin(), source_bank.end()};
#endif

    // loop over the other nodes and receive data - then write those.
    for (int i = 0; i < mpi::n_procs; ++i) {
      // number of particles for node node i
      size_t count[] {static_cast<size_t>(bank_index[i + 1] - bank_index[i])};

#ifdef OPENMC_MPI
      if (i > 0)
        MPI_Recv(source_bank.data(), count[0], mpi::source_site, i, i,
          mpi::intracomm, MPI_STATUS_IGNORE);
#endif
      // now write the source_bank data again.
      for (const auto& site : source_bank) {
        // particle is now at the iterator
        // write it to the mcpl-file
        mcpl_particle_t p;
        p.position[0] = site.r.x;
        p.position[1] = site.r.y;
        p.position[2] = site.r.z;

        // mcpl requires that the direction vector is unit length
        // which is also the case in openmc
        p.direction[0] = site.u.x;
        p.direction[1] = site.u.y;
        p.direction[2] = site.u.z;

        // MCPL stores kinetic energy in [MeV], time in [ms]
        p.ekin = site.E * 1e-6;
        p.time = site.time * 1e3;
        p.weight = site.wgt;

        switch (site.particle) {
        case ParticleType::neutron:
          p.pdgcode = 2112;
          break;
        case ParticleType::photon:
          p.pdgcode = 22;
          break;
        case ParticleType::electron:
          p.pdgcode = 11;
          break;
        case ParticleType::positron:
          p.pdgcode = -11;
          break;
        }

        mcpl_add_particle(file_id, &p);
      }
    }
#ifdef OPENMC_MPI
    // Restore state of source bank
    std::copy(temp_source.begin(), temp_source.end(), source_bank.begin());
#endif
  } else {
#ifdef OPENMC_MPI
    MPI_Send(source_bank.data(), count_size, mpi::source_site, 0, mpi::rank,
      mpi::intracomm);
#endif
  }
}
#endif

//==============================================================================

void write_mcpl_source_point(const char* filename, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index)
{
  std::string filename_(filename);
  const auto extension = get_file_extension(filename_);
  if (extension == "") {
    filename_.append(".mcpl");
  } else if (extension != "mcpl") {
    warning("write_mcpl_source_point was passed a file extension differing "
            "from .mcpl, but an mcpl file will be written.");
  }

#ifdef OPENMC_MCPL
  mcpl_outfile_t file_id;

  std::string line;
  if (mpi::master) {
    file_id = mcpl_create_outfile(filename_.c_str());
    if (VERSION_DEV) {
      line = fmt::format("OpenMC {0}.{1}.{2}-development", VERSION_MAJOR,
        VERSION_MINOR, VERSION_RELEASE);
    } else {
      line = fmt::format(
        "OpenMC {0}.{1}.{2}", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE);
    }
    mcpl_hdr_set_srcname(file_id, line.c_str());
  }

  write_mcpl_source_bank(file_id, source_bank, bank_index);

  if (mpi::master) {
    mcpl_close_outfile(file_id);
  }
#endif
}

} // namespace openmc
