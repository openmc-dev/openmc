#include "openmc/urr.h"

#include <algorithm> // any_of
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

  // Read URR tables. The HDF5 format is a little
  // different from how we want it laid out in memory.
  // This array used to be called "prob_".
  xt::xtensor<double, 3> tmp_prob;
  read_dataset(group_id, "table", tmp_prob);
  auto shape = tmp_prob.shape();

  // We separate out into two matrices (one with CDF values,
  // the other with cross section sets) in order to improve
  // contiguity of memory accesses.
  const auto n_energy = shape[0];
  const auto n_cdf_values = shape[2];
  xt::xtensor<double, 2>::shape_type new_shape = {n_energy, n_cdf_values};
  cdf_values_.resize(new_shape);
  xs_values_.resize(new_shape);

  // Now fill in the values. Using manual loops here since we might
  // not have fancy xtensor slicing code written for GPU tensors.
  // The below enum gives how URR tables are laid out in our HDF5 tables.
  enum class URRTableParam {
    CUM_PROB,
    TOTAL,
    ELASTIC,
    FISSION,
    N_GAMMA,
    HEATING
  };
  for (int i_energy = 0; i_energy < n_energy; ++i_energy) {
    for (int i_cdf = 0; i_cdf < n_cdf_values; ++i_cdf) {
      cdf_values_(i_energy, i_cdf) =
        tmp_prob(i_energy, URRTableParam::CUM_PROB, i_cdf);
      xs_values_(i_energy, i_cdf).total =
        tmp_prob(i_energy, URRTableParam::TOTAL, i_cdf);
      xs_values_(i_energy, i_cdf).elastic =
        tmp_prob(i_energy, URRTableParam::ELASTIC, i_cdf);
      xs_values_(i_energy, i_cdf).fission =
        tmp_prob(i_energy, URRTableParam::FISSION, i_cdf);
      xs_values_(i_energy, i_cdf).n_gamma =
        tmp_prob(i_energy, URRTableParam::N_GAMMA, i_cdf);
      xs_values_(i_energy, i_cdf).heating =
        tmp_prob(i_energy, URRTableParam::HEATING, i_cdf);
    }
  }
}

bool UrrData::has_negative() const
{

  // Lambda checks if any value in XSSset is negative
  auto xs_set_negative = [](const XSSet& xs) {
    return xs.total < 0.0 || xs.elastic < 0.0 || xs.fission < 0.0 ||
           xs.n_gamma < 0.0 || xs.heating < 0.0;
  };

  return std::any_of(cdf_values_.begin(), cdf_values_.end(), [](double x) {
    return x < 0.0;
  }) || std::any_of(xs_values_.begin(), xs_values_.end(), xs_set_negative);
}

} // namespace openmc
