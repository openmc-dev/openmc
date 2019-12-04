#include "openmc/mgxs_interface.h"

#include <string>
#include <unordered_set>

#include "openmc/cell.h"
#include "openmc/cross_sections.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry_aux.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"


namespace openmc {

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

namespace data {
  MgxsInterface mg;
}

MgxsInterface::MgxsInterface(const std::string& path_cross_sections,
                             const std::vector<std::string> xs_to_read,
                             const std::vector<std::vector<double>> xs_temps)
{
  read_header(path_cross_sections);
  set_nuclides_and_temperatures(xs_to_read, xs_temps);
  init();
}

void MgxsInterface::set_nuclides_and_temperatures(
    std::vector<std::string> xs_to_read,
    std::vector<std::vector<double>> xs_temps)
{
  // Check to remove all duplicates
  xs_to_read_ = xs_to_read;
  xs_temps_to_read_ = xs_temps;
  if (xs_to_read_.size() != xs_temps.size())
    fatal_error("The list of macro XS temperatures to read does not "
                "correspond in length to the number of XS names. ");
}

void MgxsInterface::init()
{

  // Check that at least some data was set to be read
  if (xs_to_read_.size() == 0)
    warning("No MGXS nuclides were set to be read.");

  // Check if MGXS Library exists
  if (!file_exists(cross_sections_path_)) {
    // Could not find MGXS Library file
    fatal_error("Cross sections HDF5 file '" + cross_sections_path_ +
      "' does not exist.");
  }

  write_message("Loading cross section data...", 5);

  // Open file for reading
  hid_t file_id = file_open(cross_sections_path_, 'r');

  // Read filetype
  std::string type;
  read_attribute(file_id, "filetype", type);
  if (type != "mgxs") {
    fatal_error("Provided MGXS Library is not a MGXS Library file.");
  }

  // Read revision number for the MGXS Library file and make sure it matches
  // with the current version
  std::array<int, 2> array;
  read_attribute(file_id, "version", array);
  if (array != VERSION_MGXS_LIBRARY) {
    fatal_error("MGXS Library file version does not match current version "
      "supported by OpenMC.");
  }

  // ==========================================================================
  // READ ALL MGXS CROSS SECTION TABLES
  for (unsigned i_nuc=0; i_nuc<xs_to_read_.size(); ++i_nuc)
    add_mgxs(file_id, xs_to_read_[i_nuc], xs_temps_to_read_[i_nuc]);

  file_close(file_id);

  create_macro_xs();
}

//==============================================================================

void
MgxsInterface::add_mgxs(hid_t file_id, const std::string& name,
  const std::vector<double>& temperature)
{
  write_message("Loading " + std::string(name) + " data...", 6);

  // Check to make sure cross section set exists in the library
  hid_t xs_grp;
  if (object_exists(file_id, name.c_str())) {
    xs_grp = open_group(file_id, name.c_str());
  } else {
    fatal_error("Data for " + std::string(name) + " does not exist in "
                + "provided MGXS Library");
  }

  nuclides_.emplace_back(xs_grp, temperature, num_energy_groups_,
      num_delayed_groups_);
  close_group(xs_grp);
}

//==============================================================================

void MgxsInterface::create_macro_xs()
{
  // Get temperatures to read for each material
  auto kTs = get_mat_kTs();

  // Force all nuclides in a material to be the same representation.
  // Therefore type(nuclides[mat->nuclide_[0]]) dictates type(macroxs).
  // At the same time, we will find the scattering type, as that will dictate
  // how we allocate the scatter object within macroxs.
  for (int i = 0; i < model::materials.size(); ++i) {
    if (kTs[i].size() > 0) {
      // Convert atom_densities to a vector
      auto& mat {model::materials[i]};
      std::vector<double> atom_densities(mat->atom_density_.begin(),
        mat->atom_density_.end());

      // Build array of pointers to nuclides's Mgxs objects needed for this
      // material
      std::vector<Mgxs*> mgxs_ptr;
      for (int i_nuclide : mat->nuclide_) {
        mgxs_ptr.push_back(&nuclides_[i_nuclide]);
      }

      macro_xs_.emplace_back(mat->name_, kTs[i], mgxs_ptr, atom_densities,
          num_energy_groups_, num_delayed_groups_);
    } else {
      // Preserve the ordering of materials by including a blank entry
      macro_xs_.emplace_back();
    }
  }
}

//==============================================================================

std::vector<std::vector<double>> MgxsInterface::get_mat_kTs()
{
  std::vector<std::vector<double>> kTs(model::materials.size());

  for (const auto& cell : model::cells) {
    // Skip non-material cells
    if (cell->fill_ != C_NONE) continue;

    for (int j = 0; j < cell->material_.size(); ++j) {
      // Skip void materials
      int i_material = cell->material_[j];
      if (i_material == MATERIAL_VOID) continue;

      // Get temperature of cell (rounding to nearest integer)
      double sqrtkT = cell->sqrtkT_.size() == 1 ?
        cell->sqrtkT_[j] : cell->sqrtkT_[0];
      double kT = sqrtkT * sqrtkT;

      // Add temperature if it hasn't already been added
      if (!contains(kTs[i_material], kT)) {
        kTs[i_material].push_back(kT);
      }
    }
  }
  return kTs;
}

//==============================================================================

void MgxsInterface::read_header(const std::string& path_cross_sections)
{
  // Save name of HDF5 file to be read to struct data
  cross_sections_path_ = path_cross_sections;

  // Check if MGXS Library exists
  if (!file_exists(cross_sections_path_)) {
    // Could not find MGXS Library file
    fatal_error("Cross sections HDF5 file '" + cross_sections_path_ +
      "' does not exist.");
  }
  write_message("Reading cross sections HDF5 file...", 5);

  // Open file for reading
  hid_t file_id = file_open(cross_sections_path_, 'r', true);

  ensure_exists(file_id, "energy_groups", true);
  read_attribute(file_id, "energy_groups", num_energy_groups_);

  if (attribute_exists(file_id, "delayed_groups")) {
    read_attribute(file_id, "delayed_groups", num_delayed_groups_);
  } else {
    num_delayed_groups_ = 0;
  }

  ensure_exists(file_id, "group structure", true);
  read_attribute(file_id, "group structure", rev_energy_bins_);

  // Reverse energy bins
  std::copy(rev_energy_bins_.crbegin(), rev_energy_bins_.crend(),
    std::back_inserter(energy_bins_));

  // Create average energies
  for (int i = 0; i < energy_bins_.size() - 1; ++i) {
    energy_bin_avg_.push_back(0.5*
    (energy_bins_[i] + energy_bins_[i+1]));
  }

  // Add entries into libraries for MG data
  xs_names_ = group_names(file_id);
  if (xs_names_.empty()) {
    fatal_error("At least one MGXS data set must be present in mgxs "
      "library file!");
  }

  // Close MGXS HDF5 file
  file_close(file_id);
}

void put_mgxs_header_data_to_globals()
{
  // Get the minimum and maximum energies
  int neutron = static_cast<int>(Particle::Type::neutron);
  data::energy_min[neutron] = data::mg.energy_bins_.back();
  data::energy_max[neutron] = data::mg.energy_bins_.front();

  // Save available XS names to library list, so that when
  // materials are read, the specified mgxs can be confirmed
  // as present
  for (auto& name : data::mg.xs_names_) {
    Library lib {};
    lib.type_ = Library::Type::neutron;
    lib.materials_.push_back(name);
    data::libraries.push_back(lib);
  }
}

void set_mg_interface_nuclides_and_temps()
{
  // Get temperatures from global data
  std::vector<std::vector<double>> nuc_temps(data::nuclide_map.size());
  std::vector<std::vector<double>> dummy;
  get_temperatures(nuc_temps, dummy);

  // Build vector of nuclide names which are to be read
  std::vector<std::string> nuclide_names(data::nuclide_map.size());
  for (const auto& kv : data::nuclide_map) {
    nuclide_names[kv.second] = kv.first;
  }

  std::unordered_set<std::string> already_read;

  // Loop over materials to find xs and temperature to be read
  for (const auto& mat : model::materials) {
    for (int i_nuc : mat->nuclide_) {
      std::string& name = nuclide_names[i_nuc];

      if (already_read.find(name) == already_read.end()) {
        data::mg.xs_to_read_.push_back(name);
        data::mg.xs_temps_to_read_.push_back(nuc_temps[i_nuc]);
        already_read.insert(name);
      }
    }
  }
}

void mark_fissionable_mgxs_materials()
{
  // Loop over all files
  for (const auto& mat : model::materials) {
    for (int i_nuc : mat->nuclide_) {
      if (data::mg.nuclides_[i_nuc].fissionable) {
        mat->fissionable_ = true;
      }
    }
  }
}

} // namespace openmc
