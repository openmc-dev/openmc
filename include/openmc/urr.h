//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! UrrData contains probability tables for the unresolved resonance range.
//==============================================================================

class UrrData {
public:
  // Since we access all of these at once, we want
  // them contiguous in memory.
  struct XSSet {
    double total;
    double elastic;
    double fission;
    double n_gamma;
    double heating;
  };

  Interpolation interp_; //!< interpolation type
  int inelastic_flag_;   //!< inelastic competition flag
  int absorption_flag_;  //!< other absorption flag
  bool multiply_smooth_; //!< multiply by smooth cross section?

  vector<double> energy_; //!< incident energies
  auto n_energy() const { return energy_.size(); }

  /* The row indexes correspond to the incident energy table, and column
   * indices correspond to values of the CDF at that energy. For the CDF matrix
   * below, obviously, values of the CDF are stored. For the xs_values
   * variable, the columns line up with the index of cdf_values.
   */
  xt::xtensor<double, 2> cdf_values_; // Note: must be row major!
  xt::xtensor<XSSet, 2> xs_values_;

  // Number of points in the CDF
  auto n_cdf() const { return cdf_values_.shape()[1]; }

  //! \brief Load the URR data from the provided HDF5 group
  explicit UrrData(hid_t group_id);

  // Checks if any negative CDF or XS values are present
  bool has_negative() const;

  // Checks if the passed energy is within the bounds of the URR table
  bool energy_in_bounds(double E) const
  {
    return energy_.front() < E && E < energy_.back();
  }
};

} // namespace openmc

#endif // OPENMC_URR_H
