#include "openmc/urr.h"

namespace openmc {

void
UrrData::from_hdf5(hid_t group_id)
{
  // Read interpolation and other flags
  read_attribute(group_id, "interpolation", interp_);
  read_attribute(group_id, "inelastic", inelastic_flag_);
  read_attribute(group_id, "absorption", absorption_flag_);
  int i;
  read_attribute(group_id, "multiply_smooth", i);
  multiply_smooth_ = (i == 1);

  // read the enrgies at which tables exist
  hid_t dset = open_dataset(group_id, "energy");
  hsize_t dims[1];
  get_shape(dset, dims);
  close_dataset(dset);
  n_energy_ = static_cast<int>(dims[0]);
  energy_ = xt::xtensor<double, 1>({n_energy_}, 0.);
  read_dataset_as_shape(group_id, "energy", energy_);

  // Read URR tables
  dset = open_dataset(group_id, "table");
  hsize_t dims3[3];
  get_shape(dset, dims3);
  close_dataset(dset);
  n_prob_ = static_cast<int>(dims3[0]);
  xt::xarray<double> temp({n_energy_, 6, n_prob_});
  read_dataset(group_id, "table", temp);

  prob_ = xt::xtensor<double, 3>({n_energy_, 6, n_prob_}, 0.);
  prob_ = temp;
  //TODO: swap 1st and last indices?

}

}