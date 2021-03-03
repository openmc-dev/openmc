//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"

namespace openmc {

//==============================================================================
//! UrrData contains probability tables for the unresolved resonance range.
//==============================================================================

class UrrData{
public:
  Interpolation interp_;          //!< interpolation type
  int inelastic_flag_;            //!< inelastic competition flag
  int absorption_flag_;           //!< other absorption flag
  bool multiply_smooth_;          //!< multiply by smooth cross section?
  int n_energy_;                  //!< number of energy points
  xt::xtensor<double, 1> energy_; //!< incident energies
  xt::xtensor<double, 3> prob_;   //!< Actual probability tables

  //! \brief Load the URR data from the provided HDF5 group
  explicit UrrData(hid_t group_id);
};

} // namespace openmc

#endif // OPENMC_URR_H
