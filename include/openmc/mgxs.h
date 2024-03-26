//! \file mgxs.h
//! A collection of classes for Multi-Group Cross Section data

#ifndef OPENMC_MGXS_H
#define OPENMC_MGXS_H

#include <string>

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/vector.h"
#include "openmc/xsdata.h"

namespace openmc {

//==============================================================================
// MGXS contains the mgxs data for a nuclide/material
//==============================================================================

class Mgxs {
private:
  xt::xtensor<double, 1> kTs; // temperature in eV (k * T)
  AngleDistributionType
    scatter_format;       // flag for if this is legendre, histogram, or tabular
  int num_groups;         // number of energy groups
  int num_delayed_groups; // number of delayed neutron groups
  vector<XsData> xs;      // Cross section data
  // MGXS Incoming Flux Angular grid information
  int n_pol;
  int n_azi;
  vector<double> polar;
  vector<double> azimuthal;

  //! \brief Initializes the Mgxs object metadata
  //!
  //! @param in_name Name of the object.
  //! @param in_awr atomic-weight ratio.
  //! @param in_kTs temperatures (in units of eV) that data is available.
  //! @param in_fissionable Is this item fissionable or not.
  //! @param in_scatter_format Denotes whether Legendre, Tabular, or
  //!   Histogram scattering is used.
  //! @param in_is_isotropic Is this an isotropic or angular with respect to
  //!   the incoming particle.
  //! @param in_polar Polar angle grid.
  //! @param in_azimuthal Azimuthal angle grid.
  void init(const std::string& in_name, double in_awr,
    const vector<double>& in_kTs, bool in_fissionable,
    AngleDistributionType in_scatter_format, bool in_is_isotropic,
    const vector<double>& in_polar, const vector<double>& in_azimuthal);

  //! \brief Initializes the Mgxs object metadata from the HDF5 file
  //!
  //! @param xs_id HDF5 group id for the cross section data.
  //! @param temperature Temperatures to read.
  //! @param temps_to_read Resultant list of temperatures in the library
  //!   to read which correspond to the requested temperatures.
  //! @param order_dim Resultant dimensionality of the scattering order.
  void metadata_from_hdf5(hid_t xs_id, const vector<double>& temperature,
    vector<int>& temps_to_read, int& order_dim);

  //! \brief Performs the actual act of combining the microscopic data for a
  //!   single temperature.
  //!
  //! @param micros Microscopic objects to combine.
  //! @param scalars Scalars to multiply the microscopic data by.
  //! @param micro_ts The temperature index of the microscopic objects that
  //!   corresponds to the temperature of interest.
  //! @param this_t The temperature index of the macroscopic object.
  void combine(const vector<Mgxs*>& micros, const vector<double>& scalars,
    const vector<int>& micro_ts, int this_t);

  //! \brief Checks to see if this and that are able to be combined
  //!
  //! This comparison is used when building macroscopic cross sections
  //! from microscopic cross sections.
  //! @param that The other Mgxs to compare to this one.
  //! @return True if they can be combined, False otherwise.
  bool equiv(const Mgxs& that);

public:
  std::string name; // name of dataset, e.g., UO2
  double awr;       // atomic weight ratio
  bool fissionable; // Is this fissionable
  bool is_isotropic {
    true}; // used to skip search for angle indices if isotropic

  Mgxs() = default;

  //! \brief Constructor that loads the Mgxs object from the HDF5 file
  //!
  //! @param xs_id HDF5 group id for the cross section data.
  //! @param temperature Temperatures to read.
  //! @param num_group number of energy groups
  //! @param num_delay number of delayed groups
  Mgxs(hid_t xs_id, const vector<double>& temperature, int num_group,
    int num_delay);

  //! \brief Constructor that initializes and populates all data to build a
  //!   macroscopic cross section from microscopic cross sections.
  //!
  //! @param in_name Name of the object.
  //! @param mat_kTs temperatures (in units of eV) that data is needed.
  //! @param micros Microscopic objects to combine.
  //! @param atom_densities Atom densities of those microscopic quantities.
  //! @param num_group number of energy groups
  //! @param num_delay number of delayed groups
  Mgxs(const std::string& in_name, const vector<double>& mat_kTs,
    const vector<Mgxs*>& micros, const vector<double>& atom_densities,
    int num_group, int num_delay);

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
  //! @param t Temperature index.
  //! @param a Angle index.
  //! @return Requested cross section value.
  double get_xs(MgxsType xstype, int gin, const int* gout, const double* mu,
    const int* dg, int t, int a);

  inline double get_xs(MgxsType xstype, int gin, int t, int a)
  {
    return get_xs(xstype, gin, nullptr, nullptr, nullptr, t, a);
  }

  //! \brief Samples the fission neutron energy and if prompt or delayed.
  //!
  //! @param gin Incoming energy group.
  //! @param dg Sampled delayed group index.
  //! @param gout Sampled outgoing energy group.
  //! @param seed Pseudorandom seed pointer
  //! @param t Temperature index.
  //! @param a Angle index.
  void sample_fission_energy(
    int gin, int& dg, int& gout, uint64_t* seed, int t, int a);

  //! \brief Samples the outgoing energy and angle from a scatter event.
  //!
  //! @param gin Incoming energy group.
  //! @param gout Sampled outgoing energy group.
  //! @param mu Sampled cosine of the change-in-angle.
  //! @param wgt Weight of the particle to be adjusted.
  //! @param seed Pseudorandom seed pointer.
  //! @param t Temperature index.
  //! @param a Angle index.
  void sample_scatter(
    int gin, int& gout, double& mu, double& wgt, uint64_t* seed, int t, int a);

  //! \brief Calculates cross section quantities needed for tracking.
  //!
  //! @param p The particle whose attributes set which MGXS to get.
  void calculate_xs(Particle& p);

  //! \brief Sets the temperature index in the particle's cache.
  //!
  //! @param p Particle.
  void set_temperature_index(Particle& p);

  //! \brief Gets the temperature index given a temperature.
  //!
  //! @param sqrtkT Temperature of the material.
  //! @return The temperature index corresponding to sqrtkT.
  int get_temperature_index(double sqrtkT) const;

  //! \brief Sets the angle index in the particle's cache.
  //!
  //! @param p Particle.
  void set_angle_index(Particle& p);

  //! \brief Gets the angle index given a direction.
  //!
  //! @param u Incoming particle direction.
  //! @return The angle index corresponding to u.
  int get_angle_index(const Direction& u) const;

  //! \brief Provide const access to list of XsData held by this
  const vector<XsData>& get_xsdata() const { return xs; }
};

} // namespace openmc
#endif // OPENMC_MGXS_H
