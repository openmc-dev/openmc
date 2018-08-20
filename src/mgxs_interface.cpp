#include "openmc/mgxs_interface.h"

#include <string>

#include "openmc/error.h"
#include "openmc/math_functions.h"


namespace openmc {

//==============================================================================
// Mgxs data loading interface methods
//==============================================================================

void
add_mgxs_c(hid_t file_id, const char* name, int energy_groups,
     int delayed_groups, int n_temps, const double temps[],  double tolerance,
     int max_order, bool legendre_to_tabular, int legendre_to_tabular_points,
     int& method)
{
  // Convert temps to a vector for the from_hdf5 function
  double_1dvec temperature(temps, temps + n_temps);

  write_message("Loading " + std::string(name) + " data...", 6);

  // Check to make sure cross section set exists in the library
  hid_t xs_grp;
  if (object_exists(file_id, name)) {
    xs_grp = open_group(file_id, name);
  } else {
    fatal_error("Data for " + std::string(name) + " does not exist in "
                + "provided MGXS Library");
  }

  Mgxs mg(xs_grp, energy_groups, delayed_groups, temperature, tolerance,
          max_order, legendre_to_tabular, legendre_to_tabular_points, method);

  nuclides_MG.push_back(mg);
  close_group(xs_grp);
}

//==============================================================================

bool
query_fissionable_c(int n_nuclides, const int i_nuclides[])
{
  bool result = false;
  for (int n = 0; n < n_nuclides; n++) {
    if (nuclides_MG[i_nuclides[n] - 1].fissionable) result = true;
  }
  return result;
}

//==============================================================================

void
create_macro_xs_c(const char* mat_name, int n_nuclides, const int i_nuclides[],
     int n_temps, const double temps[], const double atom_densities[],
     double tolerance, int& method)
{
  if (n_temps > 0) {
    // // Convert temps to a vector
    double_1dvec temperature(temps, temps + n_temps);

    // Convert atom_densities to a vector
    double_1dvec atom_densities_vec(atom_densities,
         atom_densities + n_nuclides);

    // Build array of pointers to nuclides_MG's Mgxs objects needed for this
    // material
    std::vector<Mgxs*> mgxs_ptr(n_nuclides);
    for (int n = 0; n < n_nuclides; n++) {
      mgxs_ptr[n] = &nuclides_MG[i_nuclides[n] - 1];
    }

    Mgxs macro(mat_name, temperature, mgxs_ptr, atom_densities_vec,
         tolerance, method);
    macro_xs.emplace_back(macro);
  } else {
    // Preserve the ordering of materials by including a blank entry
    Mgxs macro;
    macro_xs.emplace_back(macro);
  }
}

//==============================================================================
// Mgxs tracking/transport/tallying interface methods
//==============================================================================

void
calculate_xs_c(int i_mat, int gin, double sqrtkT, const double uvw[3],
     double& total_xs, double& abs_xs, double& nu_fiss_xs)
{
  macro_xs[i_mat - 1].calculate_xs(gin - 1, sqrtkT, uvw, total_xs, abs_xs,
       nu_fiss_xs);
}

//==============================================================================

void
sample_scatter_c(int i_mat, int gin, int& gout, double& mu, double& wgt,
                 double uvw[3])
{
  int gout_c = gout - 1;
  macro_xs[i_mat - 1].sample_scatter(gin - 1, gout_c, mu, wgt);

  // adjust return value for fortran indexing
  gout = gout_c + 1;

  // Rotate the angle
  rotate_angle_c(uvw, mu, nullptr);
}

//==============================================================================

void
sample_fission_energy_c(int i_mat, int gin, int& dg, int& gout)
{
  int dg_c = 0;
  int gout_c = 0;
  macro_xs[i_mat - 1].sample_fission_energy(gin - 1, dg_c, gout_c);

  // adjust return values for fortran indexing
  dg = dg_c + 1;
  gout = gout_c + 1;
}

//==============================================================================

double
get_nuclide_xs_c(int index, int xstype, int gin, int* gout, double* mu, int* dg)
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
  return nuclides_MG[index - 1].get_xs(xstype, gin - 1, gout_c_p, mu, dg_c_p);
}

//==============================================================================

double
get_macro_xs_c(int index, int xstype, int gin, int* gout, double* mu, int* dg)
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
  return macro_xs[index - 1].get_xs(xstype, gin - 1, gout_c_p, mu, dg_c_p);
}

//==============================================================================

void
set_nuclide_angle_index_c(int index, const double uvw[3])
{
  // Update the values
  nuclides_MG[index - 1].set_angle_index(uvw);
}

//==============================================================================

void
set_macro_angle_index_c(int index, const double uvw[3])
{
  // Update the values
  macro_xs[index - 1].set_angle_index(uvw);
}

//==============================================================================

void
set_nuclide_temperature_index_c(int index, double sqrtkT)
{
  // Update the values
  nuclides_MG[index - 1].set_temperature_index(sqrtkT);
}

//==============================================================================
// General Mgxs methods
//==============================================================================

void
get_name_c(int index, int name_len, char* name)
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
get_awr_c(int index)
{
  return nuclides_MG[index - 1].awr;
}

} // namespace openmc
