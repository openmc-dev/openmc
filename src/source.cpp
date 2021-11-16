#include "openmc/source.h"

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#define HAS_DYNAMIC_LINKING
#endif

#include <algorithm> // for move

#ifdef HAS_DYNAMIC_LINKING
#include <dlfcn.h> // for dlopen, dlsym, dlclose, dlerror
#endif

#include <fmt/core.h>
#include "xtensor/xadapt.hpp"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
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

#ifdef __CUDACC__
#include "openmc/cuda/calculate_xs.h" // gpu::materials
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

vector<unique_ptr<Source>> external_sources;

}

#ifdef __CUDACC__
namespace gpu {
__constant__ unique_ptr<Source>* external_sources;
__constant__ int n_external_sources;
} // namespace gpu
#endif

//==============================================================================
// IndependentSource implementation
//==============================================================================

IndependentSource::IndependentSource(
  UPtrSpace&& space, UPtrAngle&& angle, unique_ptr<Distribution>&& energy)
  : space_ {std::move(space)}, angle_ {std::move(angle)}, energy_ {
                                                            std::move(energy)}
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
        space_ = UPtrSpace {make_unique<CartesianIndependent>(node_space)};
      } else if (type == "cylindrical") {
        space_ = UPtrSpace {make_unique<CylindricalIndependent>(node_space)};
      } else if (type == "spherical") {
        space_ = UPtrSpace {make_unique<SphericalIndependent>(node_space)};
      } else if (type == "box") {
        space_ = UPtrSpace {make_unique<SpatialBox>(node_space)};
      } else if (type == "fission") {
        space_ = UPtrSpace {make_unique<SpatialBox>(node_space, true)};
      } else if (type == "point") {
        space_ = UPtrSpace {make_unique<SpatialPoint>(node_space)};
      } else {
        fatal_error(fmt::format(
          "Invalid spatial distribution for external source: {}", type));
      }

    } else {
      // If no spatial distribution specified, make it a point source
      space_ = UPtrSpace {make_unique<SpatialPoint>()};
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
        angle_ = make_unique<Isotropic>();
      } else if (type == "monodirectional") {
        angle_ = make_unique<Monodirectional>(node_angle);
      } else if (type == "mu-phi") {
        angle_ = make_unique<PolarAzimuthal>(node_angle);
      } else {
        fatal_error(fmt::format(
          "Invalid angular distribution for external source: {}", type));
      }

    } else {
      angle_ = make_unique<Isotropic>();
    }

    // Determine external source energy distribution
    if (check_for_node(node, "energy")) {
      pugi::xml_node node_dist = node.child("energy");
      energy_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
      energy_ = make_unique<Watt>(0.988e6, 2.249e-6);
    }
  }
}

SourceSite IndependentSource::sample(uint64_t* seed, Particle* p) const
{
#ifdef __CUDA_ARCH__
  using gpu::cells;
  using gpu::materials;
#else
  using model::cells;
  using model::materials;
#endif

  SourceSite site;

  // Set weight to one by default
  site.wgt = 1.0;

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
    double xyz[3] {site.r.x, site.r.y, site.r.z};

    // If not provided with a particle, we have to allocate on-the-fly
    if (!p) {
#ifdef __CUDA_ARCH__
      __trap();
#else
      int err = openmc_find_cell(xyz, &cell_index, &instance);
      found = (err != OPENMC_E_GEOMETRY);
#endif
    } else {
      // Otherwise, we can re-use pre-allocated particle data
      Position pos(xyz);
      Direction dir(0, 0, 1); // just happens to be what openmc_find_cell does
      p->r() = pos;
      p->u() = dir;
      found = exhaustive_find_cell(*p);
      cell_index = p->coord(p->n_coord() - 1).cell;
      instance = p->cell_instance();
    }

    // Check if spatial site is in fissionable material
    if (found) {
      if (space_->only_fissionable()) {
        // Determine material
        const auto& c = cells[cell_index];
        auto mat_index =
          c->material_.size() == 1 ? c->material_[0] : c->material_[instance];

        if (mat_index == MATERIAL_VOID) {
          found = false;
        } else {
          if (!materials[mat_index]->fissionable_)
            found = false;
        }
      }
    } else { // !found
      ++n_reject;
      if (n_reject >= EXTSRC_REJECT_THRESHOLD &&
          static_cast<double>(n_accept)/n_reject <= EXTSRC_REJECT_FRACTION) {
#ifdef __CUDA_ARCH__
        printf("More than 95%% of sources rejected...\n");
        __trap();
#else
        fatal_error("More than 95% of external source sites sampled were "
                    "rejected. Please check your external source definition.");
#endif
      }
    }
  }

  // Increment number of accepted samples
  ++n_accept;

  // Sample angle
  site.u = angle_->sample(seed);

#ifndef __CUDA_ARCH__ // just going to not do this error checking on device
  auto ptype = static_cast<int>(particle_);
  // Check for monoenergetic source above maximum particle energy
  auto energy_ptr = dynamic_cast<Discrete*>(energy_.get());
  double const& min_energy = data::energy_min[ptype];
  double const& max_energy = data::energy_max[ptype];
  if (energy_ptr) {
    auto const& energies = energy_ptr->x();
    if (std::any_of(energies.begin(), energies.end(),
          [max_energy](double const& energy) { return energy > max_energy; })) {
      fatal_error("Source energy above range of energies of at least "
                  "one cross section table");
    } else if (std::any_of(energies.begin(), energies.end(),
                 [min_energy](
                   double const& energy) { return energy < min_energy; })) {
      fatal_error("Source energy below range of energies of at least "
                  "one cross section table");
    }
  }
#else
  if (particle_ == ParticleType::photon)
    printf("TODO make source sampling of photons work");
#endif

  while (true) {
    // Sample energy spectrum
    site.E = energy_->sample(seed);

    // Resample if energy falls outside minimum or maximum particle energy
#ifdef __CUDA_ARCH__ // Remove this after addressing the above TODO
    if (site.E < gpu::energy_max_neutron && site.E > gpu::energy_min_neutron)
      break;
#else
    if (site.E < data::energy_max[ptype] && site.E > data::energy_min[ptype])
      break;
#endif
  }

  // Set delayed group
  site.delayed_group = 0;
  // Set surface ID
  site.surf_id = 0;

  return site;
}

//==============================================================================
// FileSource implementation
//==============================================================================

FileSource::FileSource(std::string const& path)
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

SourceSite FileSource::sample(uint64_t* seed, Particle* p) const
{
  size_t i_site = sites_.size()*prn(seed);
  return sites_[i_site];
}

//==============================================================================
// CustomSourceWrapper implementation
//==============================================================================

#ifndef __CUDACC__
CustomSourceWrapper::CustomSourceWrapper(
  std::string const& path, std::string const& parameters)
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
    std::string error_msg = fmt::format("Couldn't open the openmc_create_source symbol: {}", dlsym_error);
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
  if (custom_source_.get()) custom_source_.reset();

#ifdef HAS_DYNAMIC_LINKING
  dlclose(shared_library_);
#else
  fatal_error("Custom source libraries have not yet been implemented for "
              "non-POSIX systems");
#endif
}
#endif

//==============================================================================
// Non-member functions
//==============================================================================

void initialize_source()
{
  write_message("Initializing source particles...", 5);

  if (simulation::work_per_rank != settings::max_particles_in_flight)
    fatal_error("TODO fix initialize_source to work with different number of "
                "particles in flight");

  // Generation source sites from specified distribution in user input
  #pragma omp parallel for
  for (int64_t i = 0; i < simulation::work_per_rank; ++i) {

    // This saves some time re-allocating memory if we're
    // doing a structure-of-array particle data layout.
    // sample_external_source is able to use the pre-allocated
    // space here.
    Particle* p = nullptr;
#ifdef __CUDACC__
    Particle this_part(i);
    p = &this_part;
#endif

    // initialize random number seed
    int64_t id = simulation::total_gen*settings::n_particles +
      simulation::work_index[mpi::rank] + i + 1;
    uint64_t seed = init_seed(id, STREAM_SOURCE);

    // sample external source distribution
    simulation::source_bank[i] = sample_external_source(&seed, p);
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

SourceSite sample_external_source(uint64_t* seed, Particle* p)
{
#ifdef __CUDA_ARCH__
  using gpu::external_sources;
  using gpu::n_external_sources;
#else
  using model::external_sources;
  auto n_external_sources = model::external_sources.size();
#endif

  // Determine total source strength
  double total_strength = 0.0;
  for (int i = 0; i < n_external_sources; ++i) {
    total_strength += external_sources[i]->strength();
  }

  // Sample from among multiple source distributions
  int i = 0;
  if (n_external_sources > 1) {
    double xi = prn(seed)*total_strength;
    double c = 0.0;
    for (; i < n_external_sources; ++i) {
      c += external_sources[i]->strength();
      if (xi < c) break;
    }
  }

  // Sample source site from i-th source distribution
  SourceSite site {external_sources[i]->sample(seed, p)};

  // If running in MG, convert site.E to group
#ifndef __CUDACC__
  if (!settings::run_CE) {
    site.E = lower_bound_index(data::mg.rev_energy_bins_.begin(),
      data::mg.rev_energy_bins_.end(), site.E);
    site.E = data::mg.num_energy_groups_ - site.E - 1.;
  }
#endif

  return site;
}

void free_memory_source()
{
  model::external_sources.clear();
}

} // namespace openmc
