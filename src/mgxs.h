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
  // last angle that corresponds to a
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

    //! \brief Initializes the Mgxs object metadata from the HDF5 file
    //!
    //! @param xs_id HDF5 group id for the cross section data.
    //! @param in_num_groups Number of energy groups.
    //! @param in_num_delayed_groups Number of delayed groups.
    //! @param temperature Temperatures to read.
    //! @param method Method of choosing nearest temperatures.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param temps_to_read Resultant list of temperatures in the library
    //!   to read which correspond to the requested temperatures.
    //! @param order_dim Resultant dimensionality of the scattering order.
    //! @param n_threads Number of threads at runtime.
    void
    metadata_from_hdf5(const hid_t xs_id, const int in_num_groups,
         const int in_num_delayed_groups, double_1dvec& temperature,
         int& method, const double tolerance, int_1dvec& temps_to_read,
         int& order_dim, const int n_threads);

    //! \brief Initializes the Mgxs object metadata
    //!
    //! @param in_name Name of the object.
    //! @param in_awr atomic-weight ratio.
    //! @param in_kTs temperatures (in units of eV) that data is available.
    //! @param in_fissionable Is this item fissionable or not.
    //! @param in_scatter_format Denotes whether Legendre, Tabular, or
    //!   Histogram scattering is used.
    //! @param in_num_groups Number of energy groups.
    //! @param in_num_delayed_groups Number of delayed groups.
    //! @param in_is_isotropic Is this an isotropic or angular with respect to
    //!   the incoming particle.
    //! @param in_polar Polar angle grid.
    //! @param in_azimuthal Azimuthal angle grid.
    //! @param n_threads Number of threads at runtime.
    void
    init(const std::string& in_name, const double in_awr,
         const double_1dvec& in_kTs, const bool in_fissionable,
         const int in_scatter_format, const int in_num_groups,
         const int in_num_delayed_groups, const bool in_is_isotropic,
         const double_1dvec& in_polar, const double_1dvec& in_azimuthal,
         const int n_threads);

    //! \brief Performs the actual act of combining the microscopic data for a
    //!   single temperature.
    //!
    //! @param micros Microscopic objects to combine.
    //! @param scalars Scalars to multiply the microscopic data by.
    //! @param micro_ts The temperature index of the microscopic objects that
    //!   corresponds to the temperature of interest.
    //! @param this_t The temperature index of the macroscopic object.
    void
    combine(std::vector<Mgxs*>& micros, double_1dvec& scalars,
         int_1dvec& micro_ts, int this_t);

    //! \brief Checks to see if this and that are able to be combined
    //!
    //! This comparison is used when building macroscopic cross sections
    //! from microscopic cross sections.
    //! @param that The other Mgxs to compare to this one.
    //! @return True if they can be combined, False otherwise.
    bool equiv(const Mgxs& that);

  public:

    std::string name;   // name of dataset, e.g., UO2
    double awr;         // atomic weight ratio
    bool fissionable;   // Is this fissionable
    std::vector<CacheData> cache; // index and data cache

    //! \brief Initializes and populates all data to build a macroscopic
    //!   cross section from microscopic cross section.
    //!
    //! @param in_name Name of the object.
    //! @param mat_kTs temperatures (in units of eV) that data is needed.
    //! @param micros Microscopic objects to combine.
    //! @param atom_densities Atom densities of those microscopic quantities.
    //! @param method Method of choosing nearest temperatures.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param n_threads Number of threads at runtime.
    void
    build_macro(const std::string& in_name, double_1dvec& mat_kTs,
         std::vector<Mgxs*>& micros, double_1dvec& atom_densities,
         int& method, const double tolerance, const int n_threads);

    //! \brief Loads the Mgxs object from the HDF5 file
    //!
    //! @param xs_id HDF5 group id for the cross section data.
    //! @param energy_groups Number of energy groups.
    //! @param delayed_groups Number of delayed groups.
    //! @param temperature Temperatures to read.
    //! @param method Method of choosing nearest temperatures.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param max_order Maximum order requested by the user;
    //!   this is only used for Legendre scattering.
    //! @param legendre_to_tabular Flag to denote if any Legendre provided
    //!   should be converted to a Tabular representation.
    //! @param legendre_to_tabular_points If a conversion is requested, this
    //!   provides the number of points to use in the tabular representation.
    //! @param n_threads Number of threads at runtime.
    void
    from_hdf5(hid_t xs_id, const int energy_groups,
         const int delayed_groups, double_1dvec& temperature, int& method,
         const double tolerance, const int max_order,
         const bool legendre_to_tabular, const int legendre_to_tabular_points,
         const int n_threads);

    //! \brief Provides a cross section value given certain parameters
    //!
    //! @param xstype Type of cross section requested, according to the
    //!   enumerated constants.
    //! @param gin Incoming energy group.
    //! @param gout Outgoing energy group; use nullptr if irrelevant, or if a
    //!   sum is requested.
    //! @param mu Cosine of the change-in-angle, for scattering quantities;
    //!   use nullptr if irrelevant.
    //! @param dg delayed group index; use nullptr if irrelevant.
    //! @return Requested cross section value.
    double
    get_xs(const int tid, const int xstype, const int gin, int* gout,
         double* mu, int* dg);

    //! \brief Samples the fission neutron energy and if prompt or delayed.
    //!
    //! @param tid Thread id to use when using the index cache.
    //! @param gin Incoming energy group.
    //! @param dg Sampled delayed group index.
    //! @param gout Sampled outgoing energy group.
    void
    sample_fission_energy(const int tid, const int gin, int& dg, int& gout);

    //! \brief Samples the outgoing energy and angle from a scatter event.
    //!
    //! @param tid Thread id to use when using the index cache.
    //! @param gin Incoming energy group.
    //! @param gout Sampled outgoing energy group.
    //! @param mu Sampled cosine of the change-in-angle.
    //! @param wgt Weight of the particle to be adjusted.
    void
    sample_scatter(const int tid, const int gin, int& gout, double& mu,
         double& wgt);

    //! \brief Calculates cross section quantities needed for tracking.
    //!
    //! @param tid Thread id to use when using the index cache.
    //! @param gin Incoming energy group.
    //! @param sqrtkT Temperature of the material.
    //! @param uvw Incoming particle direction.
    //! @param total_xs Resultant total cross section.
    //! @param abs_xs Resultant absorption cross section.
    //! @param nu_fiss_xs Resultant nu-fission cross section.
    void
    calculate_xs(const int tid, const int gin, const double sqrtkT,
         const double uvw[3], double& total_xs, double& abs_xs,
         double& nu_fiss_xs);

    //! \brief Sets the temperature index in cache given a temperature
    //!
    //! @param tid Thread id to use when setting the index cache.
    //! @param sqrtkT Temperature of the material.
    void
    set_temperature_index(const int tid, const double sqrtkT);

    //! \brief Sets the angle index in cache given a direction
    //!
    //! @param tid Thread id to use when setting the index cache.
    //! @param uvw Incoming particle direction.
    void
    set_angle_index(const int tid, const double uvw[3]);
};

} // namespace openmc
#endif // MGXS_H