#include "openmc/urr.h"

#include <iostream>

namespace openmc {

void
UrrData::from_hdf5(hid_t group_id)
{
  // Read interpolation and other flags
  int interp_temp;
  read_attribute(group_id, "interpolation", interp_temp);
  switch (interp_temp) {
  case static_cast<int>(Interpolation::histogram):
    interp_ = Interpolation::histogram;
    break;
  case static_cast<int>(Interpolation::lin_lin):
    interp_ = Interpolation::lin_lin;
    break;
  case static_cast<int>(Interpolation::lin_log):
    interp_ = Interpolation::lin_log;
    break;
  case static_cast<int>(Interpolation::log_lin):
    interp_ = Interpolation::log_lin;
    break;
  case static_cast<int>(Interpolation::log_log):
    interp_ = Interpolation::log_log;
  }

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