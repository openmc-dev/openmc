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

struct MgxsInterface
{
  int num_energy_groups;
  int num_delayed_groups;

  std::vector<Mgxs> nuclides_MG;
  std::vector<Mgxs> macro_xs;

  std::vector<double> energy_bins;
  std::vector<double> energy_bin_avg;
  std::vector<double> rev_energy_bins;

  MgxsInterface() = default;

  // Construct from path to cross sections file
  MgxsInterface(const std::string& path_cross_sections);

  void init(const std::string& path_cross_sections);

  void add_mgxs(hid_t file_id, const std::string& name,
         const std::vector<double>& temperature);

  void create_macro_xs();

  std::vector<std::vector<double>> get_mat_kTs();

  void read_mg_cross_sections_header();
};

namespace data {
  extern MgxsInterface mgInterface;
}

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

extern "C" void
calculate_xs_c(int i_mat, int gin, double sqrtkT, Direction u,
     double& total_xs, double& abs_xs, double& nu_fiss_xs);

double
get_nuclide_xs(int index, int xstype, int gin, const int* gout,
  const double* mu, const int* dg);

inline double
get_nuclide_xs(int index, int xstype, int gin)
{return get_nuclide_xs(index, xstype, gin, nullptr, nullptr, nullptr);}

double
get_macro_xs(int index, int xstype, int gin, const int* gout,
  const double* mu, const int* dg);

inline double
get_macro_xs(int index, int xstype, int gin)
{return get_macro_xs(index, xstype, gin, nullptr, nullptr, nullptr);}

//==============================================================================
// General Mgxs methods
//==============================================================================

extern "C" void
get_name_c(int index, int name_len, char* name);

extern "C" double
get_awr_c(int index);

} // namespace openmc
#endif // OPENMC_MGXS_INTERFACE_H
