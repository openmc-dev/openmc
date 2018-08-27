//! \file mgxs.h
//! A collection of classes for Multi-Group Cross Section data

#ifndef OPENMC_MGXS_H
#define OPENMC_MGXS_H

#include <string>
#include <vector>

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/xsdata.h"


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
    void
    init(const std::string& in_name, double in_awr, const double_1dvec& in_kTs,
         bool in_fissionable, int in_scatter_format, int in_num_groups,
         int in_num_delayed_groups, bool in_is_isotropic,
         const double_1dvec& in_polar, const double_1dvec& in_azimuthal);

    //! \brief Initializes the Mgxs object metadata from the HDF5 file
    //!
    //! @param xs_id HDF5 group id for the cross section data.
    //! @param in_num_groups Number of energy groups.
    //! @param in_num_delayed_groups Number of delayed groups.
    //! @param temperature Temperatures to read.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param temps_to_read Resultant list of temperatures in the library
    //!   to read which correspond to the requested temperatures.
    //! @param order_dim Resultant dimensionality of the scattering order.
    //! @param method Method of choosing nearest temperatures.
    void
    metadata_from_hdf5(hid_t xs_id, int in_num_groups,
         int in_num_delayed_groups, const double_1dvec& temperature,
         double tolerance, int_1dvec& temps_to_read, int& order_dim,
         int& method);

    //! \brief Performs the actual act of combining the microscopic data for a
    //!   single temperature.
    //!
    //! @param micros Microscopic objects to combine.
    //! @param scalars Scalars to multiply the microscopic data by.
    //! @param micro_ts The temperature index of the microscopic objects that
    //!   corresponds to the temperature of interest.
    //! @param this_t The temperature index of the macroscopic object.
    void
    combine(const std::vector<Mgxs*>& micros, const double_1dvec& scalars,
         const int_1dvec& micro_ts, int this_t);

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

    Mgxs() = default;

    //! \brief Constructor that loads the Mgxs object from the HDF5 file
    //!
    //! @param xs_id HDF5 group id for the cross section data.
    //! @param energy_groups Number of energy groups.
    //! @param delayed_groups Number of delayed groups.
    //! @param temperature Temperatures to read.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param max_order Maximum order requested by the user;
    //!   this is only used for Legendre scattering.
    //! @param legendre_to_tabular Flag to denote if any Legendre provided
    //!   should be converted to a Tabular representation.
    //! @param legendre_to_tabular_points If a conversion is requested, this
    //!   provides the number of points to use in the tabular representation.
    //! @param method Method of choosing nearest temperatures.
    Mgxs(hid_t xs_id, int energy_groups,
         int delayed_groups, const double_1dvec& temperature, double tolerance,
         int max_order, bool legendre_to_tabular,
         int legendre_to_tabular_points, int& method);

    //! \brief Constructor that initializes and populates all data to build a
    //!   macroscopic cross section from microscopic cross section.
    //!
    //! @param in_name Name of the object.
    //! @param mat_kTs temperatures (in units of eV) that data is needed.
    //! @param micros Microscopic objects to combine.
    //! @param atom_densities Atom densities of those microscopic quantities.
    //! @param tolerance Tolerance of temperature selection method.
    //! @param method Method of choosing nearest temperatures.
    Mgxs(const std::string& in_name, const double_1dvec& mat_kTs,
         const std::vector<Mgxs*>& micros, const double_1dvec& atom_densities,
         double tolerance, int& method);

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
    get_xs(int xstype, int gin, int* gout, double* mu, int* dg);

    //! \brief Samples the fission neutron energy and if prompt or delayed.
    //!
    //! @param gin Incoming energy group.
    //! @param dg Sampled delayed group index.
    //! @param gout Sampled outgoing energy group.
    void
    sample_fission_energy(int gin, int& dg, int& gout);

    //! \brief Samples the outgoing energy and angle from a scatter event.
    //!
    //! @param gin Incoming energy group.
    //! @param gout Sampled outgoing energy group.
    //! @param mu Sampled cosine of the change-in-angle.
    //! @param wgt Weight of the particle to be adjusted.
    void
    sample_scatter(int gin, int& gout, double& mu, double& wgt);

    //! \brief Calculates cross section quantities needed for tracking.
    //!
    //! @param gin Incoming energy group.
    //! @param sqrtkT Temperature of the material.
    //! @param uvw Incoming particle direction.
    //! @param total_xs Resultant total cross section.
    //! @param abs_xs Resultant absorption cross section.
    //! @param nu_fiss_xs Resultant nu-fission cross section.
    void
    calculate_xs(int gin, double sqrtkT, const double uvw[3],
         double& total_xs, double& abs_xs, double& nu_fiss_xs);

    //! \brief Sets the temperature index in cache given a temperature
    //!
    //! @param sqrtkT Temperature of the material.
    void
    set_temperature_index(double sqrtkT);

    //! \brief Sets the angle index in cache given a direction
    //!
    //! @param uvw Incoming particle direction.
    void
    set_angle_index(const double uvw[3]);
};

} // namespace openmc
#endif // OPENMC_MGXS_H
