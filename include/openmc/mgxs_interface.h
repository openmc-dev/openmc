//! \file mgxs_interface.h
//! A collection of C interfaces to the C++ Mgxs class

#ifndef OPENMC_MGXS_INTERFACE_H
#define OPENMC_MGXS_INTERFACE_H

#include "hdf5_interface.h"
#include "mgxs.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

extern std::vector<Mgxs> nuclides_MG;
extern std::vector<Mgxs> macro_xs;

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

extern "C" void
add_mgxs_c(hid_t file_id, const char* name, int energy_groups,
     int delayed_groups, int n_temps, const double temps[], double tolerance,
     int max_order, bool legendre_to_tabular, int legendre_to_tabular_points,
     int& method);

extern "C" bool
query_fissionable_c(int n_nuclides, const int i_nuclides[]);

extern "C" void
create_macro_xs_c(const char* mat_name, int n_nuclides, const int i_nuclides[],
     int n_temps, const double temps[], const double atom_densities[],
     double tolerance, int& method);

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

extern "C" void
calculate_xs_c(int i_mat, int gin, double sqrtkT, const double uvw[3],
     double& total_xs, double& abs_xs, double& nu_fiss_xs);

extern "C" void
sample_scatter_c(int i_mat, int gin, int& gout, double& mu, double& wgt,
     double uvw[3]);

extern "C" void
sample_fission_energy_c(int i_mat, int gin, int& dg, int& gout);

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
