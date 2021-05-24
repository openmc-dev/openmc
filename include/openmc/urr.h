//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include <iostream>
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
  double* device_energy_;
  xt::xtensor<double, 3> prob_;   //!< Actual probability tables
  double* device_prob_;
  int n_bands_;
  int n_total_prob_;

  //! \brief Load the URR data from the provided HDF5 group
  explicit UrrData(hid_t group_id);

  void flatten_urr_data();
  #pragma omp declare target
  int get_offset(int i_energy, URRTableParam type) const;
  #pragma omp end declare target
};

} // namespace openmc

#endif // OPENMC_URR_H
