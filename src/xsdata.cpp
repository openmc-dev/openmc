#include "openmc/xsdata.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>

#include "xtensor/xview.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xbuilder.hpp"

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/mgxs_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"


namespace openmc {

//==============================================================================
// XsData class methods
//==============================================================================

XsData::XsData(bool fissionable, AngleDistributionType scatter_format, int n_pol, int n_azi,
               size_t n_groups, size_t n_d_groups) :
  n_g_(n_groups),
  n_dg_(n_d_groups)
{
  size_t n_ang = n_pol * n_azi;

  // check to make sure scatter format is OK before we allocate
  if (scatter_format != AngleDistributionType::HISTOGRAM && scatter_format != AngleDistributionType::TABULAR &&
      scatter_format != AngleDistributionType::LEGENDRE) {
    fatal_error("Invalid scatter_format!");
  }
  // allocate all [temperature][angle][in group] quantities
  std::vector<size_t> shape {n_ang, n_g_};
  total = xt::zeros<double>(shape);
  absorption = xt::zeros<double>(shape);
  inverse_velocity = xt::zeros<double>(shape);
  if (fissionable) {
    fission = xt::zeros<double>(shape);
    nu_fission = xt::zeros<double>(shape);
    prompt_nu_fission = xt::zeros<double>(shape);
    kappa_fission = xt::zeros<double>(shape);
  }

  // allocate decay_rate; [temperature][angle][delayed group]
  shape[1] = n_dg_;
  decay_rate = xt::zeros<double>(shape);

  if (fissionable) {
    shape = {n_ang, n_dg_, n_g_};
    // allocate delayed_nu_fission; [temperature][angle][delay group][in group]
    delayed_nu_fission = xt::zeros<double>(shape);

    // chi_prompt; [temperature][angle][in group][out group]
    shape = {n_ang, n_g_, n_g_};
    chi_prompt = xt::zeros<double>(shape);

    // chi_delayed; [temperature][angle][delay group][in group][out group]
    shape = {n_ang, n_dg_, n_g_, n_g_};
    chi_delayed = xt::zeros<double>(shape);
  }


  for (int a = 0; a < n_ang; a++) {
    if (scatter_format == AngleDistributionType::HISTOGRAM) {
      scatter.emplace_back(new ScattDataHistogram);
    } else if (scatter_format == AngleDistributionType::TABULAR) {
      scatter.emplace_back(new ScattDataTabular);
    } else if (scatter_format == AngleDistributionType::LEGENDRE) {
      scatter.emplace_back(new ScattDataLegendre);
    }
  }
}

//==============================================================================

void
XsData::from_hdf5(hid_t xsdata_grp, bool fissionable, AngleDistributionType scatter_format,
  AngleDistributionType final_scatter_format, int order_data, bool is_isotropic, int n_pol, int n_azi)
{
  // Reconstruct the dimension information so it doesn't need to be passed
  size_t n_ang = n_pol * n_azi;
  size_t energy_groups = total.shape()[1];

  // Set the fissionable-specific data
  if (fissionable) {
    fission_from_hdf5(xsdata_grp, n_ang, is_isotropic);
  }
  // Get the non-fission-specific data
  read_nd_vector(xsdata_grp, "decay rate", decay_rate);
  read_nd_vector(xsdata_grp, "absorption", absorption, true);
  read_nd_vector(xsdata_grp, "inverse-velocity", inverse_velocity);

  // Get scattering data
  scatter_from_hdf5(xsdata_grp, n_ang, scatter_format,
    final_scatter_format, order_data);

  // Check absorption to ensure it is not 0 since it is often the
  // denominator in tally methods
  xt::filtration(absorption, xt::equal(absorption, 0.)) = 1.e-10;

  // Get or calculate the total x/s
  if (object_exists(xsdata_grp, "total")) {
    read_nd_vector(xsdata_grp, "total", total);
  } else {
    for (size_t a = 0; a < n_ang; a++) {
      for (size_t gin = 0; gin < energy_groups; gin++) {
        total(a, gin) = absorption(a, gin) + scatter[a]->scattxs[gin];
      }
    }
  }

  // Fix if total is 0, since it is in the denominator when tallying
  xt::filtration(total, xt::equal(total, 0.)) = 1.e-10;
}

//==============================================================================

void
XsData::fission_vector_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang,
  bool is_isotropic)
{
  // Data is provided as nu-fission and chi with a beta for delayed info

  // Get chi
  xt::xtensor<double, 2> temp_chi({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  temp_chi /= xt::view(xt::sum(temp_chi, {1}), xt::all(), xt::newaxis());

  // Now every incoming group in prompt_chi and delayed_chi is the normalized
  // chi we just made
  chi_prompt = xt::view(temp_chi, xt::all(), xt::newaxis(), xt::all());
  chi_delayed = xt::view(temp_chi, xt::all(), xt::newaxis(), xt::newaxis(),
                         xt::all());

  // Get nu-fission
  xt::xtensor<double, 2> temp_nufiss({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "nu-fission", temp_nufiss, true);

  // Get beta (strategy will depend upon the number of dimensions in beta)
  hid_t beta_dset = open_dataset(xsdata_grp, "beta");
  int beta_ndims = dataset_ndims(beta_dset);
  close_dataset(beta_dset);
  int ndim_target = 1;
  if (!is_isotropic) ndim_target += 2;
  if (beta_ndims == ndim_target) {
    xt::xtensor<double, 2> temp_beta({n_ang, n_dg_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // Set prompt_nu_fission = (1. - beta_total)*nu_fission
    prompt_nu_fission = temp_nufiss * (1. - xt::sum(temp_beta, {1}));

    // Set delayed_nu_fission as beta * nu_fission
    delayed_nu_fission =
         xt::view(temp_beta, xt::all(), xt::all(), xt::newaxis()) *
         xt::view(temp_nufiss, xt::all(), xt::newaxis(), xt::all());
  } else if (beta_ndims == ndim_target + 1) {
    xt::xtensor<double, 3> temp_beta({n_ang, n_dg_, n_g_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // Set prompt_nu_fission = (1. - beta_total)*nu_fission
    prompt_nu_fission = temp_nufiss * (1. - xt::sum(temp_beta, {1}));

    // Set delayed_nu_fission as beta * nu_fission
    delayed_nu_fission = temp_beta *
         xt::view(temp_nufiss, xt::all(), xt::newaxis(), xt::all());
  }
}

void
XsData::fission_vector_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get chi-prompt
  xt::xtensor<double, 2> temp_chi_p({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi-prompt", temp_chi_p, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  temp_chi_p /= xt::view(xt::sum(temp_chi_p, {1}), xt::all(), xt::newaxis());

  // Get chi-delayed
  xt::xtensor<double, 3> temp_chi_d({n_ang, n_dg_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi-delayed", temp_chi_d, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  temp_chi_d /= xt::view(xt::sum(temp_chi_d, {2}),
       xt::all(), xt::all(), xt::newaxis());

  // Now assign the prompt and delayed chis by replicating for each incoming group
  chi_prompt = xt::view(temp_chi_p, xt::all(), xt::newaxis(), xt::all());
  chi_delayed = xt::view(temp_chi_d, xt::all(), xt::all(), xt::newaxis(),
                         xt::all());

  // Get prompt and delayed nu-fission directly
  read_nd_vector(xsdata_grp, "prompt-nu-fission", prompt_nu_fission,
       true);
  read_nd_vector(xsdata_grp, "delayed-nu-fission",
       delayed_nu_fission, true);
}

void
XsData::fission_vector_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get chi
  xt::xtensor<double, 2> temp_chi({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  temp_chi /= xt::view(xt::sum(temp_chi, {1}), xt::all(), xt::newaxis());

  // Now every incoming group in self.chi is the normalized chi we just made
  chi_prompt = xt::view(temp_chi, xt::all(), xt::newaxis(), xt::all());

  // Get nu-fission directly
  read_nd_vector(xsdata_grp, "nu-fission", prompt_nu_fission, true);
}

//==============================================================================

void
XsData::fission_matrix_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Data is provided as nu-fission and chi with a beta for delayed info

  // Get nu-fission matrix
  xt::xtensor<double, 3> temp_matrix({n_ang, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "nu-fission", temp_matrix, true);

  // Get beta (strategy will depend upon the number of dimensions in beta)
  hid_t beta_dset = open_dataset(xsdata_grp, "beta");
  int beta_ndims = dataset_ndims(beta_dset);
  close_dataset(beta_dset);
  int ndim_target = 1;
  if (!is_isotropic) ndim_target += 2;
  if (beta_ndims == ndim_target) {
    xt::xtensor<double, 2> temp_beta({n_ang, n_dg_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    xt::xtensor<double, 1> temp_beta_sum({n_ang}, 0.);
    temp_beta_sum = xt::sum(temp_beta, {1});

    // prompt_nu_fission is the sum of this matrix over outgoing groups and
    // multiplied by (1 - beta_sum)
    prompt_nu_fission = xt::sum(temp_matrix, {2}) * (1. - temp_beta_sum);

    // Store chi-prompt
    chi_prompt = xt::view(1.0 - temp_beta_sum, xt::all(), xt::newaxis(),
                          xt::newaxis()) * temp_matrix;

    // delayed_nu_fission is the sum of this matrix over outgoing groups and
    // multiplied by beta
    delayed_nu_fission =
         xt::view(temp_beta, xt::all(), xt::all(), xt::newaxis()) *
         xt::view(xt::sum(temp_matrix, {2}), xt::all(), xt::newaxis(), xt::all());

    // Store chi-delayed
    chi_delayed =
         xt::view(temp_beta, xt::all(), xt::all(), xt::newaxis(), xt::newaxis()) *
         xt::view(temp_matrix, xt::all(), xt::newaxis(), xt::all(), xt::all());

  } else if (beta_ndims == ndim_target + 1) {
    xt::xtensor<double, 3> temp_beta({n_ang, n_dg_, n_g_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    xt::xtensor<double, 2> temp_beta_sum({n_ang, n_g_}, 0.);
    temp_beta_sum = xt::sum(temp_beta, {1});

    // prompt_nu_fission is the sum of this matrix over outgoing groups and
    // multiplied by (1 - beta_sum)
    prompt_nu_fission = xt::sum(temp_matrix, {2}) * (1. - temp_beta_sum);

    // Store chi-prompt
    chi_prompt = xt::view(1.0 - temp_beta_sum, xt::all(), xt::all(),
                          xt::newaxis()) * temp_matrix;

    // delayed_nu_fission is the sum of this matrix over outgoing groups and
    // multiplied by beta
    delayed_nu_fission = temp_beta *
         xt::view(xt::sum(temp_matrix, {2}), xt::all(), xt::newaxis(), xt::all());

    // Store chi-delayed
    chi_delayed =
         xt::view(temp_beta, xt::all(), xt::all(), xt::all(), xt::newaxis()) *
         xt::view(temp_matrix, xt::all(), xt::newaxis(), xt::all(), xt::all());
  }

  //Normalize both chis
  chi_prompt /= xt::view(xt::sum(chi_prompt, {2}),
       xt::all(), xt::all(), xt::newaxis());

  chi_delayed /= xt::view(xt::sum(chi_delayed, {3}),
       xt::all(), xt::all(), xt::all(), xt::newaxis());
}

void
XsData::fission_matrix_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get the prompt nu-fission matrix
  xt::xtensor<double, 3> temp_matrix_p({n_ang, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "prompt-nu-fission", temp_matrix_p, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = xt::sum(temp_matrix_p, {2});

  // chi_prompt is this matrix but normalized over outgoing groups, which we
  // have already stored in prompt_nu_fission
  chi_prompt = temp_matrix_p /
       xt::view(prompt_nu_fission, xt::all(), xt::all(), xt::newaxis());

  // Get the delayed nu-fission matrix
  xt::xtensor<double, 4> temp_matrix_d({n_ang, n_dg_, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_matrix_d, true);

  // delayed_nu_fission is the sum over outgoing groups
  delayed_nu_fission = xt::sum(temp_matrix_d, {3});

  // chi_prompt is this matrix but normalized over outgoing groups, which we
  // have already stored in prompt_nu_fission
  chi_delayed = temp_matrix_d /
       xt::view(delayed_nu_fission, xt::all(), xt::all(), xt::all(), xt::newaxis());
}

void
XsData::fission_matrix_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get nu-fission matrix
  xt::xtensor<double, 3> temp_matrix({n_ang, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "nu-fission", temp_matrix, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = xt::sum(temp_matrix, {2});

  // chi_prompt is this matrix but normalized over outgoing groups, which we
  // have already stored in prompt_nu_fission
  chi_prompt = temp_matrix / xt::view(prompt_nu_fission, xt::all(), xt::all(),
                                      xt::newaxis());
}

//==============================================================================

void
XsData::fission_from_hdf5(hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Get the fission and kappa_fission data xs; these are optional
  read_nd_vector(xsdata_grp, "fission", fission);
  read_nd_vector(xsdata_grp, "kappa-fission", kappa_fission);

  // Get the data; the strategy for doing so depends on if the data is provided
  // as a nu-fission matrix or a set of chi and nu-fission vectors
  if (object_exists(xsdata_grp, "chi") ||
      object_exists(xsdata_grp, "chi-prompt")) {
    if (n_dg_ == 0) {
      fission_vector_no_delayed_from_hdf5(xsdata_grp, n_ang);
    } else {
      if (object_exists(xsdata_grp, "beta")) {
        fission_vector_beta_from_hdf5(xsdata_grp, n_ang, is_isotropic);
      } else {
        fission_vector_no_beta_from_hdf5(xsdata_grp, n_ang);
      }
    }
  } else {
    if (n_dg_ == 0) {
      fission_matrix_no_delayed_from_hdf5(xsdata_grp, n_ang);
    } else {
      if (object_exists(xsdata_grp, "beta")) {
        fission_matrix_beta_from_hdf5(xsdata_grp, n_ang, is_isotropic);
      } else {
        fission_matrix_no_beta_from_hdf5(xsdata_grp, n_ang);
      }
    }
  }

  // Combine prompt_nu_fission and delayed_nu_fission into nu_fission
  if (n_dg_ == 0) {
    nu_fission = prompt_nu_fission;
  } else {
    nu_fission = prompt_nu_fission + xt::sum(delayed_nu_fission, {1});
  }
}

//==============================================================================

void
XsData::scatter_from_hdf5(hid_t xsdata_grp, size_t n_ang, AngleDistributionType scatter_format,
    AngleDistributionType final_scatter_format, int order_data)
{
  if (!object_exists(xsdata_grp, "scatter_data")) {
    fatal_error("Must provide scatter_data group!");
  }
  hid_t scatt_grp = open_group(xsdata_grp, "scatter_data");

  // Get the outgoing group boundary indices
  xt::xtensor<int, 2> gmin({n_ang, n_g_}, 0.);
  read_nd_vector(scatt_grp, "g_min", gmin, true);
  xt::xtensor<int, 2> gmax({n_ang, n_g_}, 0.);
  read_nd_vector(scatt_grp, "g_max", gmax, true);

  // Make gmin and gmax start from 0 vice 1 as they do in the library
  gmin -= 1;
  gmax -= 1;

  // Now use this info to find the length of a vector to hold the flattened
  // data.
  size_t length = order_data * xt::sum(gmax - gmin + 1)();

  double_4dvec input_scatt(n_ang, double_3dvec(n_g_));
  xt::xtensor<double, 1> temp_arr({length}, 0.);
  read_nd_vector(scatt_grp, "scatter_matrix", temp_arr, true);

  // Compare the number of orders given with the max order of the problem;
  // strip off the superfluous orders if needed
  int order_dim;
  if (scatter_format == AngleDistributionType::LEGENDRE) {
    order_dim = std::min(order_data - 1, settings::max_order) + 1;
  } else {
    order_dim = order_data;
  }

  // convert the flattened temp_arr to a jagged array for passing to
  // scatt data
  size_t temp_idx = 0;
  for (size_t a = 0; a < n_ang; a++) {
    for (size_t gin = 0; gin < n_g_; gin++) {
      input_scatt[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
      for (size_t i_gout = 0; i_gout < input_scatt[a][gin].size(); i_gout++) {
        input_scatt[a][gin][i_gout].resize(order_dim);
        for (size_t l = 0; l < order_dim; l++) {
          input_scatt[a][gin][i_gout][l] = temp_arr[temp_idx++];
        }
        // Adjust index for the orders we didnt take
        temp_idx += (order_data - order_dim);
      }
    }
  }

  // Get multiplication matrix
  double_3dvec temp_mult(n_ang, double_2dvec(n_g_));
  if (object_exists(scatt_grp, "multiplicity_matrix")) {
    temp_arr.resize({length / order_data});
    read_nd_vector(scatt_grp, "multiplicity_matrix", temp_arr);

    // convert the flat temp_arr to a jagged array for passing to scatt data
    size_t temp_idx = 0;
    for (size_t a = 0; a < n_ang; a++) {
      for (size_t gin = 0; gin < n_g_; gin++) {
        temp_mult[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
        for (size_t i_gout = 0; i_gout < temp_mult[a][gin].size(); i_gout++) {
          temp_mult[a][gin][i_gout] = temp_arr[temp_idx++];
        }
      }
    }
  } else {
    // Use a default: multiplicities are 1.0.
    for (size_t a = 0; a < n_ang; a++) {
      for (size_t gin = 0; gin < n_g_; gin++) {
        temp_mult[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
        for (size_t i_gout = 0; i_gout < temp_mult[a][gin].size(); i_gout++) {
          temp_mult[a][gin][i_gout] = 1.;
        }
      }
    }
  }
  close_group(scatt_grp);

  // Finally, convert the Legendre data to tabular, if needed
  if (scatter_format == AngleDistributionType::LEGENDRE &&
      final_scatter_format == AngleDistributionType::TABULAR) {
    for (size_t a = 0; a < n_ang; a++) {
      ScattDataLegendre legendre_scatt;
      xt::xtensor<int, 1> in_gmin = xt::view(gmin, a, xt::all());
      xt::xtensor<int, 1> in_gmax = xt::view(gmax, a, xt::all());

      legendre_scatt.init(in_gmin, in_gmax,
                          temp_mult[a], input_scatt[a]);

      // Now create a tabular version of legendre_scatt
      convert_legendre_to_tabular(legendre_scatt,
           *static_cast<ScattDataTabular*>(scatter[a].get()));

      scatter_format = final_scatter_format;
    }
  } else {
    // We are sticking with the current representation
    // Initialize the ScattData object with this data
    for (size_t a = 0; a < n_ang; a++) {
      xt::xtensor<int, 1> in_gmin = xt::view(gmin, a, xt::all());
      xt::xtensor<int, 1> in_gmax = xt::view(gmax, a, xt::all());
      scatter[a]->init(in_gmin, in_gmax, temp_mult[a], input_scatt[a]);
    }
  }
}

//==============================================================================

void
XsData::combine(const std::vector<XsData*>& those_xs,
                const std::vector<double>& scalars)
{
  // Combine the non-scattering data
  for (size_t i = 0; i < those_xs.size(); i++) {
    XsData* that = those_xs[i];
    if (!equiv(*that)) fatal_error("Cannot combine the XsData objects!");
    double scalar = scalars[i];
    total += scalar * that->total;
    absorption += scalar * that->absorption;
    if (i == 0) {
      inverse_velocity = that->inverse_velocity;
    }
    if (that->prompt_nu_fission.shape()[0] > 0) {
      nu_fission += scalar * that->nu_fission;
      prompt_nu_fission += scalar * that->prompt_nu_fission;
      kappa_fission += scalar * that->kappa_fission;
      fission += scalar * that->fission;
      delayed_nu_fission += scalar * that->delayed_nu_fission;
      chi_prompt += scalar * that->chi_prompt;
      chi_delayed += scalar * that->chi_delayed;
    }
    decay_rate += scalar * that->decay_rate;
  }

  // Allow the ScattData object to combine itself
  for (size_t a = 0; a < total.shape()[0]; a++) {
    // Build vector of the scattering objects to incorporate
    std::vector<ScattData*> those_scatts(those_xs.size());
    for (size_t i = 0; i < those_xs.size(); i++) {
      those_scatts[i] = those_xs[i]->scatter[a].get();
    }

    // Now combine these guys
    scatter[a]->combine(those_scatts, scalars);
  }
}

//==============================================================================

bool
XsData::equiv(const XsData& that)
{
  return (absorption.shape() == that.absorption.shape());
}

} //namespace openmc
