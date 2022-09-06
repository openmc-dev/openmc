#include "openmc/source.h"

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define HAS_DYNAMIC_LINKING
#endif

#include <algorithm> // for move

#ifdef HAS_DYNAMIC_LINKING
#include <dlfcn.h> // for dlopen, dlsym, dlclose, dlerror
#endif

#include "xtensor/xadapt.hpp"
#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/xml_interface.h"

#ifdef OPENMC_MCPL
#include <mcpl.h>
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

#ifdef OPENMC_MCPL
const bool MCPL_ENABLED = true;
#else
const bool MCPL_ENABLED = false;
#endif

namespace model {

vector<unique_ptr<Source>> external_sources;
}

//==============================================================================
// IndependentSource implementation
//==============================================================================

IndependentSource::IndependentSource(
  UPtrSpace space, UPtrAngle angle, UPtrDist energy, UPtrDist time)
  : space_ {std::move(space)}, angle_ {std::move(angle)},
    energy_ {std::move(energy)}, time_ {std::move(time)}
{}

IndependentSource::IndependentSource(pugi::xml_node node)
{
  // Check for particle type
  if (check_for_node(node, "particle")) {
    auto temp_str = get_node_value(node, "particle", true, true);
    if (temp_str == "neutron") {
      particle_ = ParticleType::neutron;
    } else if (temp_str == "photon") {
      particle_ = ParticleType::photon;
      settings::photon_transport = true;
    } else {
      fatal_error(std::string("Unknown source particle type: ") + temp_str);
    }
  }

  // Check for source strength
  if (check_for_node(node, "strength")) {
    strength_ = std::stod(get_node_value(node, "strength"));
  }

  // Check for external source file
  if (check_for_node(node, "file")) {

  } else {

    // Spatial distribution for external source
    if (check_for_node(node, "space")) {
      // Get pointer to spatial distribution
      pugi::xml_node node_space = node.child("space");

      // Check for type of spatial distribution and read
      std::string type;
      if (check_for_node(node_space, "type"))
        type = get_node_value(node_space, "type", true, true);
      if (type == "cartesian") {
        space_ = UPtrSpace {new CartesianIndependent(node_space)};
      } else if (type == "cylindrical") {
        space_ = UPtrSpace {new CylindricalIndependent(node_space)};
      } else if (type == "spherical") {
        space_ = UPtrSpace {new SphericalIndependent(node_space)};
      } else if (type == "box") {
        space_ = UPtrSpace {new SpatialBox(node_space)};
      } else if (type == "fission") {
        space_ = UPtrSpace {new SpatialBox(node_space, true)};
      } else if (type == "point") {
        space_ = UPtrSpace {new SpatialPoint(node_space)};
      } else {
        fatal_error(fmt::format(
          "Invalid spatial distribution for external source: {}", type));
      }

    } else {
      // If no spatial distribution specified, make it a point source
      space_ = UPtrSpace {new SpatialPoint()};
    }

    // Determine external source angular distribution
    if (check_for_node(node, "angle")) {
      // Get pointer to angular distribution
      pugi::xml_node node_angle = node.child("angle");

      // Check for type of angular distribution
      std::string type;
      if (check_for_node(node_angle, "type"))
        type = get_node_value(node_angle, "type", true, true);
      if (type == "isotropic") {
        angle_ = UPtrAngle {new Isotropic()};
      } else if (type == "monodirectional") {
        angle_ = UPtrAngle {new Monodirectional(node_angle)};
      } else if (type == "mu-phi") {
        angle_ = UPtrAngle {new PolarAzimuthal(node_angle)};
      } else {
        fatal_error(fmt::format(
          "Invalid angular distribution for external source: {}", type));
      }

    } else {
      angle_ = UPtrAngle {new Isotropic()};
    }

    // Determine external source energy distribution
    if (check_for_node(node, "energy")) {
      pugi::xml_node node_dist = node.child("energy");
      energy_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
      energy_ = UPtrDist {new Watt(0.988e6, 2.249e-6)};
    }

    // Determine external source time distribution
    if (check_for_node(node, "time")) {
      pugi::xml_node node_dist = node.child("time");
      time_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Constant time T=0
      double T[] {0.0};
      double p[] {1.0};
      time_ = UPtrDist {new Discrete {T, p, 1}};
    }
  }
}

SourceSite IndependentSource::sample(uint64_t* seed) const
{
  SourceSite site;

  // Repeat sampling source location until a good site has been found
  bool found = false;
  int n_reject = 0;
  static int n_accept = 0;
  while (!found) {
    // Set particle type
    site.particle = particle_;

    // Sample spatial distribution
    site.r = space_->sample(seed);

    // Now search to see if location exists in geometry
    int32_t cell_index, instance;
    double xyz[] {site.r.x, site.r.y, site.r.z};
    int err = openmc_find_cell(xyz, &cell_index, &instance);
    found = (err != OPENMC_E_GEOMETRY);

    // Check if spatial site is in fissionable material
    if (found) {
      auto space_box = dynamic_cast<SpatialBox*>(space_.get());
      if (space_box) {
        if (space_box->only_fissionable()) {
          // Determine material
          const auto& c = model::cells[cell_index];
          auto mat_index =
            c->material_.size() == 1 ? c->material_[0] : c->material_[instance];

          if (mat_index == MATERIAL_VOID) {
            found = false;
          } else {
            if (!model::materials[mat_index]->fissionable_)
              found = false;
          }
        }
      }
    }

    // Check for rejection
    if (!found) {
      ++n_reject;
      if (n_reject >= EXTSRC_REJECT_THRESHOLD &&
          static_cast<double>(n_accept) / n_reject <= EXTSRC_REJECT_FRACTION) {
        fatal_error("More than 95% of external source sites sampled were "
                    "rejected. Please check your external source definition.");
      }
    }
  }

  // Increment number of accepted samples
  ++n_accept;

  // Sample angle
  site.u = angle_->sample(seed);

  // Check for monoenergetic source above maximum particle energy
  auto p = static_cast<int>(particle_);
  auto energy_ptr = dynamic_cast<Discrete*>(energy_.get());
  if (energy_ptr) {
    auto energies = xt::adapt(energy_ptr->x());
    if (xt::any(energies > data::energy_max[p])) {
      fatal_error("Source energy above range of energies of at least "
                  "one cross section table");
    } else if (xt::any(energies < data::energy_min[p])) {
      fatal_error("Source energy below range of energies of at least "
                  "one cross section table");
    }
  }

  while (true) {
    // Sample energy spectrum
    site.E = energy_->sample(seed);

    // Resample if energy falls outside minimum or maximum particle energy
    if (site.E < data::energy_max[p] && site.E > data::energy_min[p])
      break;
  }

  // Sample particle creation time
  site.time = time_->sample(seed);

  return site;
}

//==============================================================================
// FileSource implementation
//==============================================================================

FileSource::FileSource(std::string path)
{
  // Check if source file exists
  if (!file_exists(path)) {
    fatal_error(fmt::format("Source file '{}' does not exist.", path));
  }

  // Read the source from a binary file instead of sampling from some
  // assumed source distribution
  write_message(6, "Reading source file from {}...", path);

  // Open the binary file
  hid_t file_id = file_open(path, 'r', true);

  // Check to make sure this is a source file
  std::string filetype;
  read_attribute(file_id, "filetype", filetype);
  if (filetype != "source" && filetype != "statepoint") {
    fatal_error("Specified starting source file not a source file type.");
  }

  // Read in the source particles
  read_source_bank(file_id, sites_, false);

  // Close file
  file_close(file_id);
}

SourceSite FileSource::sample(uint64_t* seed) const
{
  size_t i_site = sites_.size() * prn(seed);
  return sites_[i_site];
}

//==============================================================================
// CustomSourceWrapper implementation
//==============================================================================

CustomSourceWrapper::CustomSourceWrapper(
  std::string path, std::string parameters)
{
#ifdef HAS_DYNAMIC_LINKING
  // Open the library
  shared_library_ = dlopen(path.c_str(), RTLD_LAZY);
  if (!shared_library_) {
    fatal_error("Couldn't open source library " + path);
  }

  // reset errors
  dlerror();

  // get the function to create the custom source from the library
  auto create_custom_source = reinterpret_cast<create_custom_source_t*>(
    dlsym(shared_library_, "openmc_create_source"));

  // check for any dlsym errors
  auto dlsym_error = dlerror();
  if (dlsym_error) {
    std::string error_msg = fmt::format(
      "Couldn't open the openmc_create_source symbol: {}", dlsym_error);
    dlclose(shared_library_);
    fatal_error(error_msg);
  }

  // create a pointer to an instance of the custom source
  custom_source_ = create_custom_source(parameters);

#else
  fatal_error("Custom source libraries have not yet been implemented for "
              "non-POSIX systems");
#endif
}

CustomSourceWrapper::~CustomSourceWrapper()
{
  // Make sure custom source is cleared before closing shared library
  if (custom_source_.get())
    custom_source_.reset();

#ifdef HAS_DYNAMIC_LINKING
  dlclose(shared_library_);
#else
  fatal_error("Custom source libraries have not yet been implemented for "
              "non-POSIX systems");
#endif
}

#ifdef OPENMC_MCPL
//===========================================================================
// Read particles from an MCPL-file
//===========================================================================
MCPLFileSource::MCPLFileSource(std::string path)
{
  // Check if source file exists
  if (!file_exists(path)) {
    fatal_error(fmt::format("Source file '{}' does not exist.", path));
  }

  // Read the source from a binary file instead of sampling from some
  // assumed source distribution
  write_message(6, "Reading mcpl source file from {}",path);

  // Open the mcpl file
  mcpl_file = mcpl_open_file(path.c_str());

  //do checks on the mcpl_file to see if particles are many enough.
  // should model this on the example source shown in the docs.
  n_sites=mcpl_hdr_nparticles(mcpl_file);

  read_source_bank(sites_);
}

MCPLFileSource::~MCPLFileSource(){
  mcpl_close_file(mcpl_file);
}

SourceSite MCPLFileSource::sample(uint64_t* seed) const
{
  size_t i_site = sites_.size() * prn(seed);
  return sites_[i_site];
}

void MCPLFileSource::read_source_bank(vector<SourceSite> &sites_)
{
  sites_.resize(n_sites);
  for (int i=0;i<n_sites;i++){
    sites_[i]=read_single_particle();
  }
}

SourceSite MCPLFileSource::read_single_particle() const
{
  SourceSite omc_particle_;
  const mcpl_particle_t *mcpl_particle;
  // extract particle from mcpl-file
  mcpl_particle=mcpl_read(mcpl_file);
  // check if it is a neutron, photon, electron, or positron. Otherwise skip.
  int pdg=mcpl_particle->pdgcode;
  while ( pdg!=2112 && pdg!=22 && pdg!=11 && pdg!=-11) {
    mcpl_particle=mcpl_read(mcpl_file);
    pdg=mcpl_particle->pdgcode;
    //should check for file exhaustion. This could happen if particles are other than
    //neutrons, photons, electrons, or positrons.
  }

  switch(pdg){
    case 2112:
      omc_particle_.particle=ParticleType::neutron;
      break;
    case 22:
      omc_particle_.particle=ParticleType::photon;
      break;
    case 11:
      omc_particle_.particle=ParticleType::electron;
      break;
    case -11:
      omc_particle_.particle=ParticleType::positron;
      break;
  }

  //particle is good, convert to openmc-formalism
  omc_particle_.r.x=mcpl_particle->position[0];
  omc_particle_.r.y=mcpl_particle->position[1];
  omc_particle_.r.z=mcpl_particle->position[2];

  omc_particle_.u.x=mcpl_particle->direction[0];
  omc_particle_.u.y=mcpl_particle->direction[1];
  omc_particle_.u.z=mcpl_particle->direction[2];

  //mcpl stores kinetic energy in MeV
  omc_particle_.E=mcpl_particle->ekin*1e6;
  //mcpl stores time in ms
  omc_particle_.time=mcpl_particle->time*1e-3;
  omc_particle_.wgt=mcpl_particle->weight;
  omc_particle_.delayed_group=0;
  return omc_particle_;
}
#endif //OPENMC_MCPL

//==============================================================================
// Non-member functions
//==============================================================================

void initialize_source()
{
  write_message("Initializing source particles...", 5);

// Generation source sites from specified distribution in user input
#pragma omp parallel for
  for (int64_t i = 0; i < simulation::work_per_rank; ++i) {
    // initialize random number seed
    int64_t id = simulation::total_gen * settings::n_particles +
                 simulation::work_index[mpi::rank] + i + 1;
    uint64_t seed = init_seed(id, STREAM_SOURCE);

    // sample external source distribution
    simulation::source_bank[i] = sample_external_source(&seed);
  }

  // Write out initial source
  if (settings::write_initial_source) {
    write_message("Writing out initial source...", 5);
    std::string filename = settings::path_output + "initial_source.h5";
    hid_t file_id = file_open(filename, 'w', true);
    write_source_bank(file_id, false);
    file_close(file_id);
  }
}

SourceSite sample_external_source(uint64_t* seed)
{
  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : model::external_sources)
    total_strength += s->strength();

  // Sample from among multiple source distributions
  int i = 0;
  if (model::external_sources.size() > 1) {
    double xi = prn(seed) * total_strength;
    double c = 0.0;
    for (; i < model::external_sources.size(); ++i) {
      c += model::external_sources[i]->strength();
      if (xi < c)
        break;
    }
  }

  // Sample source site from i-th source distribution
  SourceSite site {model::external_sources[i]->sample(seed)};

  // If running in MG, convert site.E to group
  if (!settings::run_CE) {
    site.E = lower_bound_index(data::mg.rev_energy_bins_.begin(),
      data::mg.rev_energy_bins_.end(), site.E);
    site.E = data::mg.num_energy_groups_ - site.E - 1.;
  }

  return site;
}

void free_memory_source()
{
  model::external_sources.clear();
}

} // namespace openmc
