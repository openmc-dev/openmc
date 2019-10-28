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
// Global variable definitions
//==============================================================================

namespace data {

int num_energy_groups;
int num_delayed_groups;
std::vector<double> energy_bins;
std::vector<double> energy_bin_avg;
std::vector<double> rev_energy_bins;

} // namesapce data

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

void read_mgxs()
{
  // Check if MGXS Library exists
  if (!file_exists(settings::path_cross_sections)) {
    // Could not find MGXS Library file
    fatal_error("Cross sections HDF5 file '" + settings::path_cross_sections +
      "' does not exist.");
  }

  write_message("Loading cross section data...", 5);

  // Get temperatures
  std::vector<std::vector<double>> nuc_temps(data::nuclide_map.size());
  std::vector<std::vector<double>> dummy;
  get_temperatures(nuc_temps, dummy);

  // Open file for reading
  hid_t file_id = file_open(settings::path_cross_sections, 'r');

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

  std::unordered_set<std::string> already_read;

  // Build vector of nuclide names
  std::vector<std::string> nuclide_names(data::nuclide_map.size());
  for (const auto& kv : data::nuclide_map) {
    nuclide_names[kv.second] = kv.first;
  }

  // Loop over all files
  for (const auto& mat : model::materials) {
    for (int i_nuc : mat->nuclide_) {
      std::string& name = nuclide_names[i_nuc];

      if (already_read.find(name) == already_read.end()) {
        add_mgxs(file_id, name, nuc_temps[i_nuc]);
        already_read.insert(name);
      }

      if (data::nuclides_MG[i_nuc].fissionable) {
        mat->fissionable_ = true;
      }
    }
  }

  file_close(file_id);
}

//==============================================================================

void
add_mgxs(hid_t file_id, const std::string& name,
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

  data::nuclides_MG.emplace_back(xs_grp, temperature);
  close_group(xs_grp);
}

//==============================================================================

void create_macro_xs()
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

      // Build array of pointers to nuclides_MG's Mgxs objects needed for this
      // material
      std::vector<Mgxs*> mgxs_ptr;
      for (int i_nuclide : mat->nuclide_) {
        mgxs_ptr.push_back(&data::nuclides_MG[i_nuclide]);
      }

      data::macro_xs.emplace_back(mat->name_, kTs[i], mgxs_ptr, atom_densities);
    } else {
      // Preserve the ordering of materials by including a blank entry
      data::macro_xs.emplace_back();
    }
  }
}

//==============================================================================

std::vector<std::vector<double>> get_mat_kTs()
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

void read_mg_cross_sections_header()
{
  // Check if MGXS Library exists
  if (!file_exists(settings::path_cross_sections)) {
    // Could not find MGXS Library file
    fatal_error("Cross sections HDF5 file '" + settings::path_cross_sections +
      "' does not exist.");
  }
  write_message("Reading cross sections HDF5 file...", 5);

  // Open file for reading
  hid_t file_id = file_open(settings::path_cross_sections, 'r', true);

  ensure_exists(file_id, "energy_groups", true);
  read_attribute(file_id, "energy_groups", data::num_energy_groups);

  if (attribute_exists(file_id, "delayed_groups")) {
    read_attribute(file_id, "delayed_groups", data::num_delayed_groups);
  } else {
    data::num_delayed_groups = 0;
  }

  ensure_exists(file_id, "group structure", true);
  read_attribute(file_id, "group structure", data::rev_energy_bins);

  // Reverse energy bins
  std::copy(data::rev_energy_bins.crbegin(), data::rev_energy_bins.crend(),
    std::back_inserter(data::energy_bins));

  // Create average energies
  for (int i = 0; i < data::energy_bins.size() - 1; ++i) {
    data::energy_bin_avg.push_back(0.5*(data::energy_bins[i] + data::energy_bins[i+1]));
  }

  // Add entries into libraries for MG data
  auto names = group_names(file_id);
  if (names.empty()) {
    fatal_error("At least one MGXS data set must be present in mgxs "
      "library file!");
  }

  for (auto& name : names) {
    Library lib {};
    lib.type_ = Library::Type::neutron;
    lib.materials_.push_back(name);
    data::libraries.push_back(lib);
  }

  // Get the minimum and maximum energies
  int neutron = static_cast<int>(Particle::Type::neutron);
  data::energy_min[neutron] = data::energy_bins.back();
  data::energy_max[neutron] = data::energy_bins.front();

  // Close MGXS HDF5 file
  file_close(file_id);
}

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

void
calculate_xs_c(int i_mat, int gin, double sqrtkT, Direction u,
     double& total_xs, double& abs_xs, double& nu_fiss_xs)
{
  data::macro_xs[i_mat].calculate_xs(gin - 1, sqrtkT, u, total_xs, abs_xs,
       nu_fiss_xs);
}

//==============================================================================

double
get_nuclide_xs(int index, int xstype, int gin, const int* gout,
  const double* mu, const int* dg)
{
  int gout_c;
  const int* gout_c_p;
  if (gout != nullptr) {
    gout_c = *gout - 1;
    gout_c_p = &gout_c;
  } else {
    gout_c_p = gout;
  }
  return data::nuclides_MG[index].get_xs(xstype, gin - 1, gout_c_p, mu, dg);
}

//==============================================================================

double
get_macro_xs(int index, int xstype, int gin, const int* gout,
  const double* mu, const int* dg)
{
  int gout_c;
  const int* gout_c_p;
  if (gout != nullptr) {
    gout_c = *gout - 1;
    gout_c_p = &gout_c;
  } else {
    gout_c_p = gout;
  }
  return data::macro_xs[index].get_xs(xstype, gin - 1, gout_c_p, mu, dg);
}

//==============================================================================
// General Mgxs methods
//==============================================================================

void
get_name_c(int index, int name_len, char* name)
{
  // First blank out our input string
  std::string str(name_len - 1, ' ');
  std::strcpy(name, str.c_str());

  // Now get the data and copy to the C-string
  str = data::nuclides_MG[index - 1].name;
  std::strcpy(name, str.c_str());

  // Finally, remove the null terminator
  name[std::strlen(name)] = ' ';
}

//==============================================================================

double
get_awr_c(int index)
{
  return data::nuclides_MG[index - 1].awr;
}

} // namespace openmc
