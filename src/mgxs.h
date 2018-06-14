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
#include <iostream>

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
    double_1dvec kTs;   // temperature in eV (k * T)
    int scatter_format; // flag for if this is legendre, histogram, or tabular
    int num_delayed_groups; // number of delayed neutron groups
    int num_groups;     // number of energy groups
    double last_sqrtkT; // cache of the temperature corresponding to index_temp
    std::vector<XsData> xs; // Cross section data
    int n_pol;
    int n_azi;
    double_1dvec polar;
    double_1dvec azimuthal;
    void _metadata_from_hdf5(const hid_t xs_id, const int in_num_groups,
         const int in_num_delayed_groups, double_1dvec& temperature,
         int& method, const double tolerance, int_1dvec& temps_to_read,
         int& order_dim, bool& is_isotropic);
    bool equiv(const Mgxs& that);

  public:
    std::string name;   // name of dataset, e.g., UO2
    double awr;         // atomic weight ratio
    bool fissionable;   // Is this fissionable
    // TODO: The following attributes be private when Fortran is fully replaced
    int index_pol; // cache for the angle indices
    int index_azi;
    double last_uvw[3]; // cache of the angle corresponding to the above indices
    int index_temp;     // cache of temperature index
    void init(const std::string& in_name, const double in_awr,
         const double_1dvec& in_kTs, const bool in_fissionable,
         const int in_scatter_format, const int in_num_groups,
         const int in_num_delayed_groups, const double_1dvec& in_polar,
         const double_1dvec& in_azimuthal);
    void build_macro(const std::string& in_name, double_1dvec& mat_kTs,
                     std::vector<Mgxs*>& micros, double_1dvec& atom_densities,
                     int& method, double tolerance);
    void combine(std::vector<Mgxs*>& micros, double_1dvec& scalars,
                 int_1dvec& micro_ts, int this_t);
    void from_hdf5(hid_t xs_id, int energy_groups, int delayed_groups,
         double_1dvec& temperature, int& method, double tolerance,
         int max_order, bool legendre_to_tabular,
         int legendre_to_tabular_points);
    double get_xs(const int xstype, const int gin, int* gout, double* mu,
                  int* dg);
    void sample_fission_energy(const int gin, int& dg, int& gout);
    void sample_scatter(const int gin, int& gout, double& mu, double& wgt);
    void calculate_xs(const int gin, const double sqrtkT, const double uvw[3],
         double& total_xs, double& abs_xs, double& nu_fiss_xs);
    void set_temperature_index(const double sqrtkT);
    void set_angle_index(const double uvw[3]);
};

} // namespace openmc
#endif // MGXS_H