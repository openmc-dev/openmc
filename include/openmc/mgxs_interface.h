//! \file mgxs_interface.h
//! A collection of C interfaces to the C++ Mgxs class

#ifndef OPENMC_MGXS_INTERFACE_H
#define OPENMC_MGXS_INTERFACE_H

#include "hdf5_interface.h"
#include "mgxs.h"

#include <vector>

namespace openmc {

//==============================================================================
// Global MGXS data container structure
//==============================================================================

class MgxsInterface {
public:

  MgxsInterface() = default;

  // Construct from path to cross sections file, as well as a list
  // of XS to read and the corresponding temperatures for each XS
  MgxsInterface(const std::string& path_cross_sections,
                const std::vector<std::string> xs_to_read,
                const std::vector<std::vector<double>> xs_temps);

  // Does things to construct after the nuclides and temperatures to
  // read have been specified.
  void init();

  // Set which nuclides and temperatures are to be read
  void set_nuclides_and_temperatures(std::vector<std::string> xs_to_read,
                                     std::vector<std::vector<double>> xs_temps);

  // Add an Mgxs object to be managed
  void add_mgxs(hid_t file_id, const std::string& name,
         const std::vector<double>& temperature);

  // Reads just the header of the cross sections file, to find
  // min & max energies as well as the available XS
  void read_header(const std::string& path_cross_sections);

  // Calculate microscopic cross sections from nuclide macro XS
  void create_macro_xs();

  // Get the kT values which are used in the OpenMC model
  std::vector<std::vector<double>> get_mat_kTs();

  int num_energy_groups_;
  int num_delayed_groups_;
  std::vector<std::string> xs_names_; // available names in HDF5 file
  std::vector<std::string> xs_to_read_; // XS which appear in materials
  std::vector<std::vector<double>> xs_temps_to_read_; // temperatures used
  std::string cross_sections_path_; // path to MGXS h5 file
  std::vector<Mgxs> nuclides_;
  std::vector<Mgxs> macro_xs_;
  std::vector<double> energy_bins_;
  std::vector<double> energy_bin_avg_;
  std::vector<double> rev_energy_bins_;
  std::vector<std::vector<double>> nuc_temps_; // all available temperatures
};

namespace data {
  extern MgxsInterface mg;
}

// Puts available XS in MGXS file to globals so that when
// materials are read, the MGXS specified in a material can
// be ensured to be present in the available data.
void put_mgxs_header_data_to_globals();

// Set which nuclides and temperatures are to be read on
// mg through global data
void set_mg_interface_nuclides_and_temps();

// After macro XS have been read, materials can be marked as fissionable
void mark_fissionable_mgxs_materials();

} // namespace openmc
#endif // OPENMC_MGXS_INTERFACE_H
