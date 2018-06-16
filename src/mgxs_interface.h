//! \file mgxs_interface.h
//! A collection of C interfaces to the C++ Mgxs class

#ifndef MGXS_INTERFACE_H
#define MGXS_INTERFACE_H

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
add_mgxs_c(hid_t file_id, char* name, const int energy_groups,
     const int delayed_groups, const int n_temps, double temps[], int& method,
     const double tolerance, const int max_order,
     const bool legendre_to_tabular, const int legendre_to_tabular_points,
     const int n_threads);

extern "C" bool
query_fissionable_c(const int n_nuclides, const int i_nuclides[]);

extern "C" void
create_macro_xs_c(char* mat_name, const int n_nuclides,
     const int i_nuclides[], const int n_temps, const double temps[],
     const double atom_densities[], int& method, const double tolerance,
     const int n_threads);

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

extern "C" void
calculate_xs_c(const int i_mat, const int tid, const int gin,
     const double sqrtkT, const double uvw[3], double& total_xs,
     double& abs_xs, double& nu_fiss_xs);

extern "C" void
sample_scatter_c(const int i_mat, const int tid, const int gin,
     int& gout, double& mu, double& wgt, double uvw[3]);

extern "C" void
sample_fission_energy_c(const int i_mat, const int tid,
     const int gin, int& dg, int& gout);

extern "C" double
get_nuclide_xs_c(const int index, const int tid,
     const int xstype, const int gin, int* gout, double* mu, int* dg);

extern "C" double
get_macro_xs_c(const int index, const int tid,
     const int xstype, const int gin, int* gout, double* mu, int* dg);

extern "C" void
set_nuclide_angle_index_c(const int index, const int tid,
     const double uvw[3]);

extern "C" void
set_macro_angle_index_c(const int index, const int tid,
     const double uvw[3]);

extern "C" void
set_nuclide_temperature_index_c(const int index, const int tid,
     const double sqrtkT);

//==============================================================================
// General Mgxs methods
//==============================================================================

extern "C" void
get_name_c(const int index, int name_len, char* name);

extern "C" double
get_awr_c(const int index);

} // namespace openmc
#endif // MGXS_INTERFACE_H