#include "openmc/mgxs_interface.h"

#include <string>

#include "openmc/error.h"
#include "openmc/math_functions.h"


namespace openmc {

//==============================================================================
// Global variable definitions
//==============================================================================

std::vector<double> energy_bins;
std::vector<double> energy_bin_avg;
std::vector<double> rev_energy_bins;

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
  std::vector<double> temperature(temps, temps + n_temps);

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
    std::vector<double> temperature(temps, temps + n_temps);

    // Convert atom_densities to a vector
    std::vector<double> atom_densities_vec(atom_densities,
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

void read_mg_cross_sections_header_c(hid_t file_id)
{
  ensure_exists(file_id, "energy_groups", true);
  read_attribute(file_id, "energy_groups", num_energy_groups);

  ensure_exists(file_id, "group structure", true);
  read_attribute(file_id, "group structure", energy_bins);

  // Create reverse energy bins
  std::copy(energy_bins.crbegin(), energy_bins.crend(),
    std::back_inserter(rev_energy_bins));

  // Create average energies
  for (int i = 0; i < energy_bins.size() - 1; ++i) {
    energy_bin_avg.push_back(0.5*(energy_bins[i] + energy_bins[i+1]));
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
  std::string str(name_len - 1, ' ');
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
