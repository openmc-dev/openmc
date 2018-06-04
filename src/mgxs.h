//! \file mgxs.h
//! A collection of classes for Multi-Group Cross Section data

#ifndef MGXS_H
#define MGXS_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <valarray>
#include <vector>

#include "constants.h"
#include "hdf5_interface.h"
#include "math_functions.h"
#include "random_lcg.h"
#include "scattdata.h"
#include "string_functions.h"
#include "xsdata.h"


namespace openmc {


//==============================================================================
// MGXS contains the mgxs data for a nuclide/material
//==============================================================================

class Mgxs {
  private:
    std::string name;   // name of dataset, e.g., UO2
    double awr;         // atomic weight ratio
    double_1dvec kTs;   // temperature in eV (k * T)
    int scatter_format; // flag for if this is legendre, histogram, or tabular
    int num_delayed_groups; // number of delayed neutron groups
    int num_groups;     // number of energy groups
    int index_temp;     // cache of temperature index
    double last_sqrtkT; // cache of the temperature corresponding to index_temp
    std::vector<XsData> xs; // Cross section data
    int n_pol;
    int n_azi;
    int index_pol; // cache fof the angle indices
    int index_azi;
    double_1dvec polar;
    double_1dvec azimuthal;
    dir_arr last_uvw;
    void _metadata_from_hdf5(hid_t xs_id, int in_num_groups,
         int in_num_delayed_groups, double_1dvec temperature, int& method,
         double tolerance, double_1dvec& temps_to_read, int& order_dim);

  public:
    bool fissionable;   // Is this fissionable
    void init(const std::string& in_name, double in_awr, double_1dvec& in_kTs,
              bool in_fissionable, int in_scatter_format, int in_num_groups,
              int in_num_delayed_groups, double_1dvec& in_polar,
              double_1dvec& in_azimuthal);
    void build_macro(const std::string& in_name, double_1dvec& mat_kTs,
                     std::vector<Mgxs>& micros, double_1dvec& atom_densities,
                     int& method, double tolerance);
    void combine(std::vector<Mgxs>& micros, double_1dvec& scalars,
                 int_1dvec& micro_ts, int this_t);
    void from_hdf5(hid_t xs_id, int energy_groups, int delayed_groups,
                           double_1dvec temperature, int& method,
                           double tolerance, int max_order,
                           bool legendre_to_tabular,
                           int legendre_to_tabular_points);
    double get_xs(const char* xstype, int gin, int* gout, double* mu,
                  int* dg);
    void sample_fission_energy(int gin, double nu_fission, int& dg, int& gout);
    void sample_scatter(dir_arr& uvw, int gin, int& gout, double& mu,
                        double& wgt);
    void calculate_xs(int gin, double sqrtkT, dir_arr& uvw, double& total_xs,
                      double& abs_xs, double& nu_fiss_xs);
    bool equiv(const Mgxs& that);
    inline void set_temperature_index(double sqrtkT);
    inline void set_angle_index(dir_arr& uvw);
};

extern "C" void read_mgxs_library(hid_t file_id, int n_nuclides, char** names,
     int energy_groups, int delayed_groups, int n_temps, double temps[],
     int& method, double tolerance, int max_order, bool legendre_to_tabular,
     int legendre_to_tabular_points);
extern "C" bool query_fissionable(const int i_nuclides[], const int n_nuclides);
void create_macro_xs(int n_materials, double_2dvec& mat_kTs,
     std::vector<std::string>& mat_names, double_1dvec& atom_densities,
     int& method, double tolerance);


// Storage for the MGXS data
std::vector<Mgxs> nuclides_MG;
std::vector<Mgxs> macro_xs;

} // namespace openmc
#endif // MGXS_H