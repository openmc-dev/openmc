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

  // read the enrgies at which tables exist
  hid_t dset = open_dataset(group_id, "energy");
  hsize_t dims[1];
  get_shape(dset, dims);
  close_dataset(dset);
  n_energy_ = static_cast<int>(dims[0]);
  energy_ = xt::xtensor<double, 1>({dims[0]}, 0.);
  read_dataset_as_shape(group_id, "energy", energy_);

  // Read URR tables
  dset = open_dataset(group_id, "table");
  hsize_t dims3[3];
  get_shape(dset, dims3);
  close_dataset(dset);
  xt::xarray<double> temp_arr({dims3[0], dims3[1], dims3[2]}, 0.);
  read_dataset(group_id, "table", temp_arr);
  prob_ = temp_arr;
}

}