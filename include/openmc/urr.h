//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "xtensor/xtensor.hpp"

#include "openmc/hdf5_interface.h"

namespace openmc {

//==============================================================================
//! UrrData contains probability tables for the unresolved resonance range.
//==============================================================================

class UrrData{
public:
  long unsigned int n_energy_;    // # of incident energies
  long unsigned int n_prob_;      // # of probabilities
  int interp_;                    // interpolation type (2=lin-lin, 5=log-log)
  int inelastic_flag_;            // inelastic competition flag
  int absorption_flag_;           // other absorption flag
  bool multiply_smooth_;          // multiply by smooth cross section?
  xt::xtensor<double, 1> energy_; // incident energies
  xt::xtensor<double, 3> prob_;   // Actual probability tables

  //! Load the URR data from the provided HDF5 group
  void
  from_hdf5(hid_t group_id);
};

} // namespace openmc

#endif // OPENMC_URR_H
