//! \file mgxs_interface.h
//! A collection of C interfaces to the C++ Mgxs class

#ifndef OPENMC_MGXS_INTERFACE_H
#define OPENMC_MGXS_INTERFACE_H

#include "hdf5_interface.h"
#include "mgxs.h"

#include <vector>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

extern std::vector<Mgxs> nuclides_MG;
extern std::vector<Mgxs> macro_xs;
extern "C" int num_energy_groups;
extern "C" int num_delayed_groups;
extern std::vector<double> energy_bins;
extern std::vector<double> energy_bin_avg;
extern std::vector<double> rev_energy_bins;

} // namespace data

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

void read_mgxs();

extern "C" void
add_mgxs_c(hid_t file_id, const std::string& name,
     const std::vector<double>& temperature);

extern "C" bool
query_fissionable_c(int n_nuclides, const int i_nuclides[]);

extern "C" void
create_macro_xs_c(const char* mat_name, int n_nuclides, const int i_nuclides[],
     int n_temps, const double temps[], const double atom_densities[],
     double tolerance, int& method);

extern "C" void read_mg_cross_sections_header_c(hid_t file_id);

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

extern "C" void
calculate_xs_c(int i_mat, int gin, double sqrtkT, const double uvw[3],
     double& total_xs, double& abs_xs, double& nu_fiss_xs);

extern "C" double
get_nuclide_xs_c(int index, int xstype, int gin, int* gout, double* mu, int* dg);

extern "C" double
get_macro_xs_c(int index, int xstype, int gin, int* gout, double* mu, int* dg);

extern "C" void
set_nuclide_angle_index_c(int index, const double uvw[3]);

extern "C" void
set_macro_angle_index_c(int index, const double uvw[3]);

extern "C" void
set_nuclide_temperature_index_c(int index, double sqrtkT);

//==============================================================================
// General Mgxs methods
//==============================================================================

extern "C" void
get_name_c(int index, int name_len, char* name);

extern "C" double
get_awr_c(int index);

} // namespace openmc
#endif // OPENMC_MGXS_INTERFACE_H
