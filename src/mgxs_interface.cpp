#include "mgxs_interface.h"

namespace openmc {

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

void
add_mgxs_c(hid_t file_id, char* name, const int energy_groups,
     const int delayed_groups, const int n_temps, double temps[], int& method,
     const double tolerance, const int max_order,
     const bool legendre_to_tabular, const int legendre_to_tabular_points,
     const int n_threads)
{
  //!! mgxs_data.F90 will be modified to just create the list of names
  //!! in the order needed
  // Convert temps to a vector for the from_hdf5 function
  double_1dvec temperature;
  temperature.assign(temps, temps + n_temps);

  // TODO: C++ replacement for write_message
  // write_message("Loading " + std::string(names[i]) + " data...", 6);

  // Check to make sure cross section set exists in the library
  hid_t xs_grp;
  if (object_exists(file_id, name)) {
    xs_grp = open_group(file_id, name);
  } else {
    fatal_error("Data for " + std::string(name) + " does not exist in "
                + "provided MGXS Library");
  }

  Mgxs mg;
  mg.from_hdf5(xs_grp, energy_groups, delayed_groups,
       temperature, method, tolerance, max_order, legendre_to_tabular,
       legendre_to_tabular_points, n_threads);

  nuclides_MG.push_back(mg);
}

//==============================================================================

bool
query_fissionable_c(const int n_nuclides, const int i_nuclides[])
{
  bool result = false;
  for (int n = 0; n < n_nuclides; n++) {
    if (nuclides_MG[i_nuclides[n] - 1].fissionable) result = true;
  }
  return result;
}

//==============================================================================

void
create_macro_xs_c(char* mat_name, const int n_nuclides,
     const int i_nuclides[], const int n_temps, const double temps[],
     const double atom_densities[], int& method, const double tolerance,
     const int n_threads)
{
  Mgxs macro;
  if (n_temps > 0) {
    // // Convert temps to a vector
    double_1dvec temperature;
    temperature.assign(temps, temps + n_temps);

    // Convert atom_densities to a vector
    double_1dvec atom_densities_vec;
    atom_densities_vec.assign(atom_densities, atom_densities + n_nuclides);

    // Build array of pointers to nuclides_MG's Mgxs objects needed for this
    // material
    std::vector<Mgxs*> mgxs_ptr(n_nuclides);
    for (int n = 0; n < n_nuclides; n++) {
      mgxs_ptr[n] = &nuclides_MG[i_nuclides[n] - 1];
    }

    macro.build_macro(mat_name, temperature, mgxs_ptr, atom_densities_vec,
         method, tolerance, n_threads);
  }
  macro_xs.push_back(macro);
}

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

void
calculate_xs_c(const int i_mat, const int tid, const int gin,
     const double sqrtkT, const double uvw[3], double& total_xs, double& abs_xs,
     double& nu_fiss_xs)
{
  macro_xs[i_mat - 1].calculate_xs(tid, gin - 1, sqrtkT, uvw, total_xs, abs_xs,
       nu_fiss_xs);
}

//==============================================================================

void
sample_scatter_c(const int i_mat, const int tid, const int gin, int& gout,
     double& mu, double& wgt, double uvw[3])
{
  int gout_c = gout - 1;
  macro_xs[i_mat - 1].sample_scatter(tid, gin - 1, gout_c, mu, wgt);

  // adjust return value for fortran indexing
  gout = gout_c + 1;

  // Rotate the angle
  rotate_angle_c(uvw, mu, nullptr);
}

//==============================================================================

void
sample_fission_energy_c(const int i_mat, const int tid, const int gin,
     int& dg, int& gout)
{
  int dg_c = 0;
  int gout_c = 0;
  macro_xs[i_mat - 1].sample_fission_energy(tid, gin - 1, dg_c, gout_c);

  // adjust return values for fortran indexing
  dg = dg_c + 1;
  gout = gout_c + 1;
}

//==============================================================================

double
get_nuclide_xs_c(const int index, const int tid, const int xstype,
     const int gin, int* gout, double* mu, int* dg)
{
  int gout_c;
  int* gout_c_p;
  int dg_c;
  int* dg_c_p;
  if (gout != nullptr) {
    gout_c = *gout - 1;
    gout_c_p = &gout_c;
  } else {
    gout_c_p = gout;
  }
  if (dg != nullptr) {
    dg_c = *dg - 1;
    dg_c_p = &dg_c;
  } else {
    dg_c_p = dg;
  }
  return nuclides_MG[index - 1].get_xs(tid, xstype, gin - 1, gout_c_p, mu,
                                       dg_c_p);
}

//==============================================================================

double
get_macro_xs_c(const int index, const int tid, const int xstype,
     const int gin, int* gout, double* mu, int* dg)
{
  int gout_c;
  int* gout_c_p;
  int dg_c;
  int* dg_c_p;
  if (gout != nullptr) {
    gout_c = *gout - 1;
    gout_c_p = &gout_c;
  } else {
    gout_c_p = gout;
  }
  if (dg != nullptr) {
    dg_c = *dg - 1;
    dg_c_p = &dg_c;
  } else {
    dg_c_p = dg;
  }
  return macro_xs[index - 1].get_xs(tid, xstype, gin - 1, gout_c_p, mu,
                                    dg_c_p);
}

//==============================================================================

void
set_nuclide_angle_index_c(const int index, const int tid,
     const double uvw[3])
{
  // Update the values
  nuclides_MG[index - 1].set_angle_index(tid, uvw);
}

//==============================================================================

void
set_macro_angle_index_c(const int index, const int tid,
     const double uvw[3])
{
  // Update the values
  macro_xs[index - 1].set_angle_index(tid, uvw);
}

//==============================================================================

void
set_nuclide_temperature_index_c(const int index, const int tid,
     const double sqrtkT)
{
  // Update the values
  nuclides_MG[index - 1].set_temperature_index(tid, sqrtkT);
}

//==============================================================================
// General Mgxs methods
//==============================================================================

void
get_name_c(const int index, int name_len, char* name)
{
  // First blank out our input string
  std::string str(name_len, ' ');
  std::strcpy(name, str.c_str());

  // Now get the data and copy to the C-string
  str = nuclides_MG[index - 1].name;
  std::strcpy(name, str.c_str());

  // Finally, remove the null terminator
  name[std::strlen(name)] = ' ';
}

//==============================================================================

double
get_awr_c(const int index)
{
  return nuclides_MG[index - 1].awr;
}

} // namespace openmc