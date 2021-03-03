#include "openmc/urr.h"

#include <iostream>

namespace openmc {

UrrData::UrrData(hid_t group_id)
{
  // Read interpolation and other flags
  int interp_temp;
  read_attribute(group_id, "interpolation", interp_temp);
  interp_ = static_cast<Interpolation>(interp_temp);

  // read the metadata
  read_attribute(group_id, "inelastic", inelastic_flag_);
  read_attribute(group_id, "absorption", absorption_flag_);
  int temp_multiply_smooth;
  read_attribute(group_id, "multiply_smooth", temp_multiply_smooth);
  multiply_smooth_ = (temp_multiply_smooth == 1);

  // read the energies at which tables exist
  read_dataset(group_id, "energy", energy_);

  // Set n_energy_
  n_energy_ = energy_.shape()[0];

  // Read URR tables
  read_dataset(group_id, "table", prob_);
}

}