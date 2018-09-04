#include "openmc/source.h"

#include <algorithm> // for move
#include <sstream> // for stringstream

#include "xtensor/xadapt.hpp"

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

std::vector<SourceDistribution> external_sources;

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
    // Copy path of source file
    settings::path_source = get_node_value(node, "file", false, true);

    // Check if source file exists
    if (!file_exists(settings::path_source)) {
      std::stringstream msg;
      msg << "Source file '" << settings::path_source <<  "' does not exist.";
      fatal_error(msg);
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
      } else if (type == "box") {
        space_ = UPtrSpace{new SpatialBox(node_space)};
      } else if (type == "fission") {
        space_ = UPtrSpace{new SpatialBox(node_space, true)};
      } else if (type == "point") {
        space_ = UPtrSpace{new SpatialPoint(node_space)};
      } else {
        std::stringstream msg;
        msg << "Invalid spatial distribution for external source: " << type;
        fatal_error(msg);
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
        std::stringstream msg;
        msg << "Invalid angular distribution for external source: " << type;
        fatal_error(msg);
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


Bank SourceDistribution::sample() const
{
  Bank site;

  // Set weight to one by default
  site.wgt = 1.0;

  // Repeat sampling source location until a good site has been found
  bool found = false;
  int n_reject = 0;
  static int n_accept = 0;
  while (!found) {
    // Set particle type
    site.particle = static_cast<int>(particle_);

    // Sample spatial distribution
    Position r = space_->sample();
    site.xyz[0] = r.x;
    site.xyz[1] = r.y;
    site.xyz[2] = r.z;

    // Now search to see if location exists in geometry
    int32_t cell_index, instance;
    int err = openmc_find_cell(site.xyz, &cell_index, &instance);
    found = (err != OPENMC_E_GEOMETRY);

    // Check if spatial site is in fissionable material
    if (found) {
      auto space_box = dynamic_cast<SpatialBox*>(space_.get());
      if (space_box) {
        if (space_box->only_fissionable()) {
          // Determine material
          auto c = cells[cell_index - 1];
          int32_t mat_index = c->material_[instance];
          auto m = materials[mat_index];

          if (mat_index == MATERIAL_VOID) {
            found = false;
          } else {
            bool fissionable;
            openmc_material_get_fissionable(mat_index + 1, &fissionable);
            if (!fissionable) found = false;
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
  Direction u = angle_->sample();
  site.uvw[0] = u.x;
  site.uvw[1] = u.y;
  site.uvw[2] = u.z;

  // Check for monoenergetic source above maximum particle energy
  auto p = static_cast<int>(particle_);
  auto energy_ptr = dynamic_cast<Discrete*>(energy_.get());
  if (energy_ptr) {
    auto energies = xt::adapt(energy_ptr->x());
    if (xt::any(energies > energy_max[p-1])) {
      fatal_error("Source energy above range of energies of at least "
                  "one cross section table");
    } else if (xt::any(energies < energy_min[p-1])) {
      fatal_error("Source energy below range of energies of at least "
                  "one cross section table");
    }
  }

  while (true) {
    // Sample energy spectrum
    site.E = energy_->sample();

    // Resample if energy falls outside minimum or maximum particle energy
    if (site.E < energy_max[p-1] && site.E > energy_min[p-1]) break;
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

  // Get pointer to source bank
  Bank* source_bank;
  int64_t n;
  openmc_source_bank(&source_bank, &n);

  if (settings::path_source != "") {
    // Read the source from a binary file instead of sampling from some
    // assumed source distribution

    std::stringstream msg;
    msg << "Reading source file from " << settings::path_source << "...";
    write_message(msg, 6);

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
    read_source_bank(file_id, work_index.data(), source_bank);

    // Close file
    file_close(file_id);

  } else {
    // Generation source sites from specified distribution in user input
    for (int64_t i = 0; i < openmc_work; ++i) {
      // initialize random number seed
      int64_t id = openmc_total_gen*settings::n_particles +
        work_index[openmc::mpi::rank] + i + 1;
      set_particle_seed(id);

      // sample external source distribution
      source_bank[i] = sample_external_source();
    }
  }

  // Write out initial source
  if (settings::write_initial_source) {
    write_message("Writing out initial source...", 5);
    std::string filename = settings::path_output + "initial_source.h5";
    hid_t file_id = file_open(filename, 'w', true);
    write_source_bank(file_id, work_index.data(), source_bank);
    file_close(file_id);
  }
}

extern "C" double* rev_energy_bins_ptr();

Bank sample_external_source()
{
  // Set the random number generator to the source stream.
  prn_set_stream(STREAM_SOURCE);

  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : external_sources)
    total_strength += s.strength();

  // Sample from among multiple source distributions
  int i = 0;
  if (external_sources.size() > 1) {
    double xi = prn()*total_strength;
    double c = 0.0;
    for (; i < external_sources.size(); ++i) {
      c += external_sources[i].strength();
      if (xi < c) break;
    }
  }

  // Sample source site from i-th source distribution
  Bank site = external_sources[i].sample();

  // If running in MG, convert site % E to group
  if (!settings::run_CE) {
    // Get pointer to rev_energy_bins array on Fortran side
    double* rev_energy_bins = rev_energy_bins_ptr();

    int n = num_energy_groups + 1;
    site.E = lower_bound_index(rev_energy_bins, rev_energy_bins + n, site.E);
    site.E = num_energy_groups - site.E;
  }

  // Set the random number generator back to the tracking stream.
  prn_set_stream(STREAM_TRACKING);

  return site;
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void free_memory_source()
{
  external_sources.clear();
}

extern "C" double total_source_strength()
{
  double strength = 0.0;
  for (const auto& s : external_sources) {
    strength += s.strength();
  }
  return strength;
}

// Needed in fill_source_bank_fixedsource
extern "C" int overall_generation();

//! Fill source bank at end of generation for fixed source simulations
extern "C" void fill_source_bank_fixedsource()
{
  if (settings::path_source.empty()) {
    // Get pointer to source bank
    Bank* source_bank;
    int64_t n;
    openmc_source_bank(&source_bank, &n);

    for (int64_t i = 0; i < openmc_work; ++i) {
      // initialize random number seed
      int64_t id = (openmc_total_gen + overall_generation()) *
        settings::n_particles + work_index[openmc::mpi::rank] + i + 1;
      set_particle_seed(id);

      // sample external source distribution
      source_bank[i] = sample_external_source();
    }
  }
}

} // namespace openmc
