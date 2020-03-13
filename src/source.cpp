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
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/capi.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

std::vector<SourceDistribution> external_sources;

}

namespace {

using sample_t = Particle::Bank (*)(uint64_t* seed);
sample_t custom_source_function;
void* custom_source_library;

}


//==============================================================================
// SourceDistribution implementation
//==============================================================================

SourceDistribution::SourceDistribution(UPtrSpace space, UPtrAngle angle, UPtrDist energy)
  : space_{std::move(space)}, angle_{std::move(angle)}, energy_{std::move(energy)} { }

SourceDistribution::SourceDistribution(pugi::xml_node node)
{
  // Check for particle type
  if (check_for_node(node, "particle")) {
    auto temp_str = get_node_value(node, "particle", true, true);
    if (temp_str == "neutron") {
      particle_ = Particle::Type::neutron;
    } else if (temp_str == "photon") {
      particle_ = Particle::Type::photon;
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
    // Copy path of source file
    settings::path_source = get_node_value(node, "file", false, true);

    // Check if source file exists
    if (!file_exists(settings::path_source)) {
      fatal_error(fmt::format("Source file '{}' does not exist.",
        settings::path_source));
    }
  } else if (check_for_node(node, "library")) {
    settings::path_source_library = get_node_value(node, "library", false, true);
    if (!file_exists(settings::path_source_library)) {
      fatal_error(fmt::format("Source library '{}' does not exist.",
        settings::path_source_library));
    }
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
        space_ = UPtrSpace{new CartesianIndependent(node_space)};
      } else if (type == "cylindrical") {
        space_ = UPtrSpace{new CylindricalIndependent(node_space)};
      } else if (type == "spherical") {
        space_ = UPtrSpace{new SphericalIndependent(node_space)};
      } else if (type == "box") {
        space_ = UPtrSpace{new SpatialBox(node_space)};
      } else if (type == "fission") {
        space_ = UPtrSpace{new SpatialBox(node_space, true)};
      } else if (type == "point") {
        space_ = UPtrSpace{new SpatialPoint(node_space)};
      } else {
        fatal_error(fmt::format(
          "Invalid spatial distribution for external source: {}", type));
      }

    } else {
      // If no spatial distribution specified, make it a point source
      space_ = UPtrSpace{new SpatialPoint()};
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
        angle_ = UPtrAngle{new Isotropic()};
      } else if (type == "monodirectional") {
        angle_ = UPtrAngle{new Monodirectional(node_angle)};
      } else if (type == "mu-phi") {
        angle_ = UPtrAngle{new PolarAzimuthal(node_angle)};
      } else {
        fatal_error(fmt::format(
          "Invalid angular distribution for external source: {}", type));
      }

    } else {
      angle_ = UPtrAngle{new Isotropic()};
    }

    // Determine external source energy distribution
    if (check_for_node(node, "energy")) {
      pugi::xml_node node_dist = node.child("energy");
      energy_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
      energy_ = UPtrDist{new Watt(0.988e6, 2.249e-6)};
    }
  }
}


Particle::Bank SourceDistribution::sample(uint64_t* seed) const
{
  Particle::Bank site;

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
    double xyz[] {site.r.x, site.r.y, site.r.z};

    // Now search to see if location exists in geometry
    int32_t cell_index, instance;
    int err = openmc_find_cell(xyz, &cell_index, &instance);
    found = (err != OPENMC_E_GEOMETRY);

    // Check if spatial site is in fissionable material
    if (found) {
      auto space_box = dynamic_cast<SpatialBox*>(space_.get());
      if (space_box) {
        if (space_box->only_fissionable()) {
          // Determine material
          const auto& c = model::cells[cell_index];
          auto mat_index = c->material_.size() == 1
            ? c->material_[0] : c->material_[instance];

          if (mat_index == MATERIAL_VOID) {
            found = false;
          } else {
            if (!model::materials[mat_index]->fissionable_) found = false;
          }
        }
      }
    }

    // Check for rejection
    if (!found) {
      ++n_reject;
      if (n_reject >= EXTSRC_REJECT_THRESHOLD &&
          static_cast<double>(n_accept)/n_reject <= EXTSRC_REJECT_FRACTION) {
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
    if (site.E < data::energy_max[p] && site.E > data::energy_min[p]) break;
  }

  // Set delayed group
  site.delayed_group = 0;

  return site;
}

//==============================================================================
// Non-member functions
//==============================================================================

void initialize_source()
{
  write_message("Initializing source particles...", 5);

  if (!settings::path_source.empty()) {
    // Read the source from a binary file instead of sampling from some
    // assumed source distribution

    write_message(fmt::format("Reading source file from {}...",
      settings::path_source), 6);

    // Open the binary file
    hid_t file_id = file_open(settings::path_source, 'r', true);

    // Read the file type
    std::string filetype;
    read_attribute(file_id, "filetype", filetype);

    // Check to make sure this is a source file
    if (filetype != "source" && filetype != "statepoint") {
      fatal_error("Specified starting source file not a source file type.");
    }

    // Read in the source bank
    read_source_bank(file_id);

    // Close file
    file_close(file_id);
  } else if (!settings::path_source_library.empty()) {

    write_message(fmt::format("Sampling from library source {}...",
      settings::path_source), 6);

    fill_source_bank_custom_source();

  } else {
    // Generation source sites from specified distribution in user input
    for (int64_t i = 0; i < simulation::work_per_rank; ++i) {
      // initialize random number seed
      int64_t id = simulation::total_gen*settings::n_particles +
        simulation::work_index[mpi::rank] + i + 1;
      uint64_t seed = init_seed(id, STREAM_SOURCE);

      // sample external source distribution
      simulation::source_bank[i] = sample_external_source(&seed);
    }
  }

  // Write out initial source
  if (settings::write_initial_source) {
    write_message("Writing out initial source...", 5);
    std::string filename = settings::path_output + "initial_source.h5";
    hid_t file_id = file_open(filename, 'w', true);
    write_source_bank(file_id);
    file_close(file_id);
  }
}

Particle::Bank sample_external_source(uint64_t* seed)
{
  // return values from custom source if using
  if (!settings::path_source_library.empty()) {
    return sample_custom_source_library(seed);
  }

  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : model::external_sources)
    total_strength += s.strength();

  // Sample from among multiple source distributions
  int i = 0;
  if (model::external_sources.size() > 1) {
    double xi = prn(seed)*total_strength;
    double c = 0.0;
    for (; i < model::external_sources.size(); ++i) {
      c += model::external_sources[i].strength();
      if (xi < c) break;
    }
  }

  // Sample source site from i-th source distribution
  Particle::Bank site {model::external_sources[i].sample(seed)};

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

void load_custom_source_library()
{
#ifdef HAS_DYNAMIC_LINKING

  // Open the library
  custom_source_library = dlopen(settings::path_source_library.c_str(), RTLD_LAZY);
  if (!custom_source_library) {
    fatal_error("Couldn't open source library " + settings::path_source_library);
  }

  // reset errors
  dlerror();

  // get the function from the library
  using sample_t = Particle::Bank (*)(uint64_t* seed);
  custom_source_function = reinterpret_cast<sample_t>(dlsym(custom_source_library, "sample_source"));

  // check for any dlsym errors
  auto dlsym_error = dlerror();
  if (dlsym_error) {
    dlclose(custom_source_library);
    fatal_error(fmt::format("Couldn't open the sample_source symbol: {}", dlsym_error));
  }

#else
  fatal_error("Custom source libraries have not yet been implemented for "
    "non-POSIX systems");
#endif
}

void close_custom_source_library()
{
  dlclose(custom_source_library);
}

Particle::Bank sample_custom_source_library(uint64_t* seed)
{
  return custom_source_function(seed);
}

void fill_source_bank_custom_source()
{
  // Load the custom library
  load_custom_source_library();

  // Generation source sites from specified distribution in the
  // library source
  for (int64_t i = 0; i < simulation::work_per_rank; ++i) {
    // initialize random number seed
    int64_t id = (simulation::total_gen + overall_generation()) *
      settings::n_particles + simulation::work_index[mpi::rank] + i + 1;
     uint64_t seed = init_seed(id, STREAM_SOURCE);

    // sample custom library source
    simulation::source_bank[i] = sample_custom_source_library(&seed);
  }

  // release the library
  close_custom_source_library();
}

} // namespace openmc
