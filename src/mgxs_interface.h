//! \file mgxs_interface.h
//! A collection of C interfaces to the C++ Mgxs class

#ifndef MGXS_INTERFACE_H
#define MGXS_INTERFACE_H

#include "mgxs.h"


namespace openmc {

extern std::vector<Mgxs> nuclides_MG;
extern std::vector<Mgxs> macro_xs;


extern "C" void add_mgxs_c(hid_t file_id, char* name, int energy_groups,
     int delayed_groups, int n_temps, double temps[], int& method,
     double tolerance, int max_order, bool legendre_to_tabular,
     int legendre_to_tabular_points);

extern "C" bool query_fissionable_c(const int n_nuclides, const int i_nuclides[]);

extern "C" void create_macro_xs_c(char* mat_name, const int n_nuclides,
     const int i_nuclides[], const int n_temps, const double temps[],
     const double atom_densities[], int& method, const double tolerance);

extern "C" void calculate_xs_c(const int i_mat, const int gin,
     const double sqrtkT, const double uvw[3], double& total_xs,
     double& abs_xs, double& nu_fiss_xs);

extern "C" void sample_scatter_c(const int i_mat, const int gin, int& gout,
     double& mu, double& wgt, double uvw[3]);

extern "C" void sample_fission_energy_c(const int i_mat, const int gin,
     int& dg, int& gout);

extern "C" void get_name_c(const int index, int name_len, char* name);

extern "C" double get_awr_c(const int index);

extern "C" double get_nuclide_xs_c(const int index, const int xstype,
     const int gin, int* gout, double* mu, int* dg);

extern "C" double get_macro_xs_c(const int index, const int xstype,
     const int gin, int* gout, double* mu, int* dg);

extern "C" void set_nuclide_angle_index_c(const int index, const double uvw[3],
     int& last_pol, int& last_azi, double last_uvw[3]);

extern "C" void reset_nuclide_angle_index_c(const int index, const int last_pol,
     const int last_azi, const double last_uvw[3]);

extern "C" void set_macro_angle_index_c(const int index, const double uvw[3],
     int& last_pol, int& last_azi, double last_uvw[3]);

extern "C" void reset_macro_angle_index_c(const int index, const int last_pol,
     const int last_azi, const double last_uvw[3]);

extern "C" int set_nuclide_temperature_index_c(const int index,
     const double sqrtkT);

extern "C" void reset_nuclide_temperature_index_c(const int index,
     const int last_temp);

} // namespace openmc
#endif // MGXS_INTERFACE_H