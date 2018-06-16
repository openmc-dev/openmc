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
// Cache contains the cached data for an MGXS object
//==============================================================================

struct CacheData {
  double sqrtkT; // last temperature corresponding to t
  int t; // temperature index
  int a; // angle index
  // last angle that corresponds to p and a
  double u;
  double v;
  double w;
};

//==============================================================================
// MGXS contains the mgxs data for a nuclide/material
//==============================================================================

class Mgxs {
  private:
    double_1dvec kTs;   // temperature in eV (k * T)
    int scatter_format; // flag for if this is legendre, histogram, or tabular
    int num_delayed_groups; // number of delayed neutron groups
    int num_groups;     // number of energy groups
    std::vector<XsData> xs; // Cross section data
    // MGXS Incoming Flux Angular grid information
    bool is_isotropic; // used to skip search for angle indices if isotropic
    int n_pol;
    int n_azi;
    double_1dvec polar;
    double_1dvec azimuthal;
    void _metadata_from_hdf5(const hid_t xs_id, const int in_num_groups,
         const int in_num_delayed_groups, double_1dvec& temperature,
         int& method, const double tolerance, int_1dvec& temps_to_read,
         int& order_dim, const int n_threads);
    bool equiv(const Mgxs& that);

  public:
    std::string name;   // name of dataset, e.g., UO2
    double awr;         // atomic weight ratio
    bool fissionable;   // Is this fissionable
    // TODO: The following attributes be private when Fortran is fully replaced
    std::vector<CacheData> cache; // index and data cache
    void init(const std::string& in_name, const double in_awr,
         const double_1dvec& in_kTs, const bool in_fissionable,
         const int in_scatter_format, const int in_num_groups,
         const int in_num_delayed_groups, const double_1dvec& in_polar,
         const double_1dvec& in_azimuthal, const int n_threads);
    void build_macro(const std::string& in_name, double_1dvec& mat_kTs,
                     std::vector<Mgxs*>& micros, double_1dvec& atom_densities,
                     int& method, const double tolerance, const int n_threads);
    void combine(std::vector<Mgxs*>& micros, double_1dvec& scalars,
                 int_1dvec& micro_ts, int this_t);
    void from_hdf5(hid_t xs_id, const int energy_groups,
         const int delayed_groups, double_1dvec& temperature, int& method,
         const double tolerance, const int max_order,
         const bool legendre_to_tabular, const int legendre_to_tabular_points,
         const int n_threads);
    double get_xs(const int tid, const int xstype, const int gin, int* gout,
         double* mu, int* dg);
    void sample_fission_energy(const int tid, const int gin, int& dg, int& gout);
    void sample_scatter(const int tid, const int gin, int& gout, double& mu,
         double& wgt);
    void calculate_xs(const int tid, const int gin, const double sqrtkT,
         const double uvw[3], double& total_xs, double& abs_xs, double& nu_fiss_xs);
    void set_temperature_index(const int tid, const double sqrtkT);
    void set_angle_index(const int tid, const double uvw[3]);
};

} // namespace openmc
#endif // MGXS_H