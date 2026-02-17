#include "openmc/xsdata.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>

#include "openmc/tensor.h"

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

XsData::XsData(bool fissionable, AngleDistributionType scatter_format,
  int n_pol, int n_azi, size_t n_groups, size_t n_d_groups)
  : n_g_(n_groups), n_dg_(n_d_groups)
{
  size_t n_ang = n_pol * n_azi;

  // check to make sure scatter format is OK before we allocate
  if (scatter_format != AngleDistributionType::HISTOGRAM &&
      scatter_format != AngleDistributionType::TABULAR &&
      scatter_format != AngleDistributionType::LEGENDRE) {
    fatal_error("Invalid scatter_format!");
  }
  // allocate all [temperature][angle][in group] quantities
  vector<size_t> shape {n_ang, n_g_};
  total = tensor::zeros<double>(shape);
  absorption = tensor::zeros<double>(shape);
  inverse_velocity = tensor::zeros<double>(shape);
  if (fissionable) {
    fission = tensor::zeros<double>(shape);
    nu_fission = tensor::zeros<double>(shape);
    prompt_nu_fission = tensor::zeros<double>(shape);
    kappa_fission = tensor::zeros<double>(shape);
  }

  // allocate decay_rate; [temperature][angle][delayed group]
  shape[1] = n_dg_;
  decay_rate = tensor::zeros<double>(shape);

  if (fissionable) {
    shape = {n_ang, n_dg_, n_g_};
    // allocate delayed_nu_fission; [temperature][angle][delay group][in group]
    delayed_nu_fission = tensor::zeros<double>(shape);

    // chi_prompt; [temperature][angle][in group][out group]
    shape = {n_ang, n_g_, n_g_};
    chi_prompt = tensor::zeros<double>(shape);

    // chi_delayed; [temperature][angle][delay group][in group][out group]
    shape = {n_ang, n_dg_, n_g_, n_g_};
    chi_delayed = tensor::zeros<double>(shape);
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

void XsData::from_hdf5(hid_t xsdata_grp, bool fissionable,
  AngleDistributionType scatter_format,
  AngleDistributionType final_scatter_format, int order_data, bool is_isotropic,
  int n_pol, int n_azi)
{
  // Reconstruct the dimension information so it doesn't need to be passed
  size_t n_ang = n_pol * n_azi;
  size_t energy_groups = total.shape(1);

  // Set the fissionable-specific data
  if (fissionable) {
    fission_from_hdf5(xsdata_grp, n_ang, is_isotropic);
  }
  // Get the non-fission-specific data
  read_nd_tensor(xsdata_grp, "decay-rate", decay_rate);
  read_nd_tensor(xsdata_grp, "absorption", absorption, true);
  read_nd_tensor(xsdata_grp, "inverse-velocity", inverse_velocity);

  // Get scattering data
  scatter_from_hdf5(
    xsdata_grp, n_ang, scatter_format, final_scatter_format, order_data);

  // Replace zero absorption values with a small number to avoid
  // division by zero in tally methods
  for (size_t i = 0; i < absorption.size(); i++)
    if (absorption.data()[i] == 0.0)
      absorption.data()[i] = 1.e-10;

  // Get or calculate the total x/s
  if (object_exists(xsdata_grp, "total")) {
    read_nd_tensor(xsdata_grp, "total", total);
  } else {
    for (size_t a = 0; a < n_ang; a++) {
      for (size_t gin = 0; gin < energy_groups; gin++) {
        total(a, gin) = absorption(a, gin) + scatter[a]->scattxs[gin];
      }
    }
  }

  // Replace zero total cross sections with a small number to avoid
  // division by zero in tally methods
  for (size_t i = 0; i < total.size(); i++)
    if (total.data()[i] == 0.0)
      total.data()[i] = 1.e-10;
}

//==============================================================================

void XsData::fission_vector_beta_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Data is provided as nu-fission and chi with a beta for delayed info

  // Get chi
  tensor::Tensor<double> temp_chi = tensor::zeros<double>({n_ang, n_g_});
  read_nd_tensor(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi so it sums to 1 over outgoing groups for each angle
  for (size_t a = 0; a < n_ang; a++) {
    tensor::View<double> row = temp_chi.slice(a);
    row /= row.sum();
  }

  // Replicate the energy spectrum across all incoming groups — the
  // spectrum is independent of the incoming neutron energy
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++)
      chi_prompt.slice(a, gin) = temp_chi.slice(a);

  // Same spectrum for delayed neutrons, replicated across delayed groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t d = 0; d < n_dg_; d++)
      for (size_t gin = 0; gin < n_g_; gin++)
        chi_delayed.slice(a, d, gin) = temp_chi.slice(a);

  // Get nu-fission
  tensor::Tensor<double> temp_nufiss = tensor::zeros<double>({n_ang, n_g_});
  read_nd_tensor(xsdata_grp, "nu-fission", temp_nufiss, true);

  // Get beta (strategy will depend upon the number of dimensions in beta)
  hid_t beta_dset = open_dataset(xsdata_grp, "beta");
  int beta_ndims = dataset_ndims(beta_dset);
  close_dataset(beta_dset);
  int ndim_target = 1;
  if (!is_isotropic)
    ndim_target += 2;
  if (beta_ndims == ndim_target) {
    tensor::Tensor<double> temp_beta = tensor::zeros<double>({n_ang, n_dg_});
    read_nd_tensor(xsdata_grp, "beta", temp_beta, true);

    // prompt_nu_fission = (1 - sum_of_beta) * nu_fission
    auto beta_sum = temp_beta.sum(1);
    for (size_t a = 0; a < n_ang; a++)
      for (size_t g = 0; g < n_g_; g++)
        prompt_nu_fission(a, g) = temp_nufiss(a, g) * (1.0 - beta_sum(a));

    // Delayed nu-fission is the outer product of the delayed neutron
    // fraction (beta) and the fission production rate (nu-fission)
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t g = 0; g < n_g_; g++)
          delayed_nu_fission(a, d, g) = temp_beta(a, d) * temp_nufiss(a, g);
  } else if (beta_ndims == ndim_target + 1) {
    tensor::Tensor<double> temp_beta =
      tensor::zeros<double>({n_ang, n_dg_, n_g_});
    read_nd_tensor(xsdata_grp, "beta", temp_beta, true);

    // prompt_nu_fission = (1 - sum_of_beta) * nu_fission
    // Here beta is energy-dependent, so sum over delayed groups (axis 1)
    auto beta_sum = temp_beta.sum(1);
    for (size_t a = 0; a < n_ang; a++)
      for (size_t g = 0; g < n_g_; g++)
        prompt_nu_fission(a, g) = temp_nufiss(a, g) * (1.0 - beta_sum(a, g));

    // Delayed nu-fission: beta is already energy-dependent [n_ang, n_dg, n_g],
    // so scale each delayed group's beta by the total nu-fission for that group
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t g = 0; g < n_g_; g++)
          delayed_nu_fission(a, d, g) = temp_beta(a, d, g) * temp_nufiss(a, g);
  }
}

void XsData::fission_vector_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get chi-prompt
  tensor::Tensor<double> temp_chi_p = tensor::zeros<double>({n_ang, n_g_});
  read_nd_tensor(xsdata_grp, "chi-prompt", temp_chi_p, true);

  // Normalize prompt chi so it sums to 1 over outgoing groups for each angle
  for (size_t a = 0; a < n_ang; a++) {
    tensor::View<double> row = temp_chi_p.slice(a);
    row /= row.sum();
  }

  // Get chi-delayed
  tensor::Tensor<double> temp_chi_d =
    tensor::zeros<double>({n_ang, n_dg_, n_g_});
  read_nd_tensor(xsdata_grp, "chi-delayed", temp_chi_d, true);

  // Normalize delayed chi so it sums to 1 over outgoing groups for each
  // angle and delayed group
  for (size_t a = 0; a < n_ang; a++)
    for (size_t d = 0; d < n_dg_; d++) {
      tensor::View<double> row = temp_chi_d.slice(a, d);
      row /= row.sum();
    }

  // Replicate the prompt spectrum across all incoming groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++)
      chi_prompt.slice(a, gin) = temp_chi_p.slice(a);

  // Replicate the delayed spectrum across all incoming groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t d = 0; d < n_dg_; d++)
      for (size_t gin = 0; gin < n_g_; gin++)
        chi_delayed.slice(a, d, gin) = temp_chi_d.slice(a, d);

  // Get prompt and delayed nu-fission directly
  read_nd_tensor(xsdata_grp, "prompt-nu-fission", prompt_nu_fission, true);
  read_nd_tensor(xsdata_grp, "delayed-nu-fission", delayed_nu_fission, true);
}

void XsData::fission_vector_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get chi
  tensor::Tensor<double> temp_chi = tensor::zeros<double>({n_ang, n_g_});
  read_nd_tensor(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi so it sums to 1 over outgoing groups for each angle
  for (size_t a = 0; a < n_ang; a++) {
    tensor::View<double> row = temp_chi.slice(a);
    row /= row.sum();
  }

  // Replicate the energy spectrum across all incoming groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++)
      chi_prompt.slice(a, gin) = temp_chi.slice(a);

  // Get nu-fission directly
  read_nd_tensor(xsdata_grp, "nu-fission", prompt_nu_fission, true);
}

//==============================================================================

void XsData::fission_matrix_beta_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Data is provided as nu-fission and chi with a beta for delayed info

  // Get nu-fission matrix
  tensor::Tensor<double> temp_matrix =
    tensor::zeros<double>({n_ang, n_g_, n_g_});
  read_nd_tensor(xsdata_grp, "nu-fission", temp_matrix, true);

  // Get beta (strategy will depend upon the number of dimensions in beta)
  hid_t beta_dset = open_dataset(xsdata_grp, "beta");
  int beta_ndims = dataset_ndims(beta_dset);
  close_dataset(beta_dset);
  int ndim_target = 1;
  if (!is_isotropic)
    ndim_target += 2;
  if (beta_ndims == ndim_target) {
    tensor::Tensor<double> temp_beta = tensor::zeros<double>({n_ang, n_dg_});
    read_nd_tensor(xsdata_grp, "beta", temp_beta, true);

    auto beta_sum = temp_beta.sum(1);
    auto matrix_gout_sum = temp_matrix.sum(2);

    // prompt_nu_fission = sum_gout(matrix) * (1 - beta_total)
    for (size_t a = 0; a < n_ang; a++)
      for (size_t g = 0; g < n_g_; g++)
        prompt_nu_fission(a, g) = matrix_gout_sum(a, g) * (1.0 - beta_sum(a));

    // chi_prompt = (1 - beta_total) * nu-fission matrix (unnormalized)
    for (size_t a = 0; a < n_ang; a++)
      for (size_t gin = 0; gin < n_g_; gin++)
        for (size_t gout = 0; gout < n_g_; gout++)
          chi_prompt(a, gin, gout) =
            (1.0 - beta_sum(a)) * temp_matrix(a, gin, gout);

    // Delayed nu-fission is the outer product of the delayed neutron
    // fraction (beta) and the total fission rate summed over outgoing groups
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t g = 0; g < n_g_; g++)
          delayed_nu_fission(a, d, g) = temp_beta(a, d) * matrix_gout_sum(a, g);

    // chi_delayed = beta * nu-fission matrix, expanded across delayed groups
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t gin = 0; gin < n_g_; gin++)
          for (size_t gout = 0; gout < n_g_; gout++)
            chi_delayed(a, d, gin, gout) =
              temp_beta(a, d) * temp_matrix(a, gin, gout);

  } else if (beta_ndims == ndim_target + 1) {
    tensor::Tensor<double> temp_beta =
      tensor::zeros<double>({n_ang, n_dg_, n_g_});
    read_nd_tensor(xsdata_grp, "beta", temp_beta, true);

    auto beta_sum = temp_beta.sum(1);
    auto matrix_gout_sum = temp_matrix.sum(2);

    // prompt_nu_fission = sum_gout(matrix) * (1 - beta_total)
    // Here beta is energy-dependent, so beta_sum is 2D [n_ang, n_g]
    for (size_t a = 0; a < n_ang; a++)
      for (size_t g = 0; g < n_g_; g++)
        prompt_nu_fission(a, g) =
          matrix_gout_sum(a, g) * (1.0 - beta_sum(a, g));

    // chi_prompt = (1 - beta_sum) * nu-fission matrix (unnormalized)
    for (size_t a = 0; a < n_ang; a++)
      for (size_t gin = 0; gin < n_g_; gin++)
        for (size_t gout = 0; gout < n_g_; gout++)
          chi_prompt(a, gin, gout) =
            (1.0 - beta_sum(a, gin)) * temp_matrix(a, gin, gout);

    // Delayed nu-fission: beta is energy-dependent [n_ang, n_dg, n_g],
    // scale by total fission rate summed over outgoing groups
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t g = 0; g < n_g_; g++)
          delayed_nu_fission(a, d, g) =
            temp_beta(a, d, g) * matrix_gout_sum(a, g);

    // chi_delayed = beta * nu-fission matrix, expanded across delayed groups
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg_; d++)
        for (size_t gin = 0; gin < n_g_; gin++)
          for (size_t gout = 0; gout < n_g_; gout++)
            chi_delayed(a, d, gin, gout) =
              temp_beta(a, d, gin) * temp_matrix(a, gin, gout);
  }

  // Normalize chi_prompt so it sums to 1 over outgoing groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++) {
      tensor::View<double> row = chi_prompt.slice(a, gin);
      row /= row.sum();
    }

  // Normalize chi_delayed so it sums to 1 over outgoing groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t d = 0; d < n_dg_; d++)
      for (size_t gin = 0; gin < n_g_; gin++) {
        tensor::View<double> row = chi_delayed.slice(a, d, gin);
        row /= row.sum();
      }
}

void XsData::fission_matrix_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get the prompt nu-fission matrix
  tensor::Tensor<double> temp_matrix_p =
    tensor::zeros<double>({n_ang, n_g_, n_g_});
  read_nd_tensor(xsdata_grp, "prompt-nu-fission", temp_matrix_p, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = temp_matrix_p.sum(2);

  // chi_prompt is the nu-fission matrix normalized over outgoing groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++)
      for (size_t gout = 0; gout < n_g_; gout++)
        chi_prompt(a, gin, gout) =
          temp_matrix_p(a, gin, gout) / prompt_nu_fission(a, gin);

  // Get the delayed nu-fission matrix
  tensor::Tensor<double> temp_matrix_d =
    tensor::zeros<double>({n_ang, n_dg_, n_g_, n_g_});
  read_nd_tensor(xsdata_grp, "delayed-nu-fission", temp_matrix_d, true);

  // delayed_nu_fission is the sum over outgoing groups
  delayed_nu_fission = temp_matrix_d.sum(3);

  // chi_delayed is the delayed nu-fission matrix normalized over outgoing
  // groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t d = 0; d < n_dg_; d++)
      for (size_t gin = 0; gin < n_g_; gin++)
        for (size_t gout = 0; gout < n_g_; gout++)
          chi_delayed(a, d, gin, gout) =
            temp_matrix_d(a, d, gin, gout) / delayed_nu_fission(a, d, gin);
}

void XsData::fission_matrix_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get nu-fission matrix
  tensor::Tensor<double> temp_matrix =
    tensor::zeros<double>({n_ang, n_g_, n_g_});
  read_nd_tensor(xsdata_grp, "nu-fission", temp_matrix, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = temp_matrix.sum(2);

  // chi_prompt is the nu-fission matrix normalized over outgoing groups
  for (size_t a = 0; a < n_ang; a++)
    for (size_t gin = 0; gin < n_g_; gin++)
      for (size_t gout = 0; gout < n_g_; gout++)
        chi_prompt(a, gin, gout) =
          temp_matrix(a, gin, gout) / prompt_nu_fission(a, gin);
}

//==============================================================================

void XsData::fission_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Get the fission and kappa_fission data xs; these are optional
  read_nd_tensor(xsdata_grp, "fission", fission);
  read_nd_tensor(xsdata_grp, "kappa-fission", kappa_fission);

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
    nu_fission = prompt_nu_fission + delayed_nu_fission.sum(1);
  }
}

//==============================================================================

void XsData::scatter_from_hdf5(hid_t xsdata_grp, size_t n_ang,
  AngleDistributionType scatter_format,
  AngleDistributionType final_scatter_format, int order_data)
{
  if (!object_exists(xsdata_grp, "scatter_data")) {
    fatal_error("Must provide scatter_data group!");
  }
  hid_t scatt_grp = open_group(xsdata_grp, "scatter_data");

  // Get the outgoing group boundary indices
  tensor::Tensor<int> gmin = tensor::zeros<int>({n_ang, n_g_});
  read_nd_tensor(scatt_grp, "g_min", gmin, true);
  tensor::Tensor<int> gmax = tensor::zeros<int>({n_ang, n_g_});
  read_nd_tensor(scatt_grp, "g_max", gmax, true);

  // Make gmin and gmax start from 0 vice 1 as they do in the library
  gmin -= 1;
  gmax -= 1;

  // Now use this info to find the length of a vector to hold the flattened
  // data.
  size_t length = order_data * (gmax - gmin + 1).sum();

  double_4dvec input_scatt(n_ang, double_3dvec(n_g_));
  tensor::Tensor<double> temp_arr = tensor::zeros<double>({length});
  read_nd_tensor(scatt_grp, "scatter_matrix", temp_arr, true);

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
    read_nd_tensor(scatt_grp, "multiplicity_matrix", temp_arr);

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
      tensor::Tensor<int> in_gmin(gmin.slice(a));
      tensor::Tensor<int> in_gmax(gmax.slice(a));

      legendre_scatt.init(in_gmin, in_gmax, temp_mult[a], input_scatt[a]);

      // Now create a tabular version of legendre_scatt
      convert_legendre_to_tabular(
        legendre_scatt, *static_cast<ScattDataTabular*>(scatter[a].get()));

      scatter_format = final_scatter_format;
    }
  } else {
    // We are sticking with the current representation
    // Initialize the ScattData object with this data
    for (size_t a = 0; a < n_ang; a++) {
      tensor::Tensor<int> in_gmin(gmin.slice(a));
      tensor::Tensor<int> in_gmax(gmax.slice(a));
      scatter[a]->init(in_gmin, in_gmax, temp_mult[a], input_scatt[a]);
    }
  }
}

//==============================================================================

void XsData::combine(
  const vector<XsData*>& those_xs, const vector<double>& scalars)
{
  // Combine the non-scattering data
  for (size_t i = 0; i < those_xs.size(); i++) {
    XsData* that = those_xs[i];
    if (!equiv(*that))
      fatal_error("Cannot combine the XsData objects!");
    double scalar = scalars[i];
    total += scalar * that->total;
    absorption += scalar * that->absorption;
    if (i == 0) {
      inverse_velocity = that->inverse_velocity;
    }
    if (!that->prompt_nu_fission.empty()) {
      nu_fission += scalar * that->nu_fission;
      prompt_nu_fission += scalar * that->prompt_nu_fission;
      kappa_fission += scalar * that->kappa_fission;
      fission += scalar * that->fission;
      delayed_nu_fission += scalar * that->delayed_nu_fission;
      // Accumulate chi_prompt weighted by total prompt nu-fission
      // (summed over energy groups) for this constituent
      {
        auto pnf_sum = that->prompt_nu_fission.sum(1);
        size_t n_ang = chi_prompt.shape(0);
        size_t n_g = chi_prompt.shape(1);
        for (size_t a = 0; a < n_ang; a++)
          for (size_t gin = 0; gin < n_g; gin++)
            for (size_t gout = 0; gout < n_g; gout++)
              chi_prompt(a, gin, gout) +=
                scalar * pnf_sum(a) * that->chi_prompt(a, gin, gout);
      }
      // Accumulate chi_delayed weighted by total delayed nu-fission
      // (summed over energy groups) for this constituent
      {
        auto dnf_sum = that->delayed_nu_fission.sum(2);
        size_t n_ang = chi_delayed.shape(0);
        size_t n_dg = chi_delayed.shape(1);
        size_t n_g = chi_delayed.shape(2);
        for (size_t a = 0; a < n_ang; a++)
          for (size_t d = 0; d < n_dg; d++)
            for (size_t gin = 0; gin < n_g; gin++)
              for (size_t gout = 0; gout < n_g; gout++)
                chi_delayed(a, d, gin, gout) +=
                  scalar * dnf_sum(a, d) * that->chi_delayed(a, d, gin, gout);
      }
    }
    decay_rate += scalar * that->decay_rate;
  }

  // Normalize chi_prompt so it sums to 1 over outgoing groups
  {
    size_t n_ang = chi_prompt.shape(0);
    size_t n_g = chi_prompt.shape(1);
    for (size_t a = 0; a < n_ang; a++)
      for (size_t gin = 0; gin < n_g; gin++) {
        tensor::View<double> row = chi_prompt.slice(a, gin);
        row /= row.sum();
      }
  }
  // Normalize chi_delayed so it sums to 1 over outgoing groups
  {
    size_t n_ang = chi_delayed.shape(0);
    size_t n_dg = chi_delayed.shape(1);
    size_t n_g = chi_delayed.shape(2);
    for (size_t a = 0; a < n_ang; a++)
      for (size_t d = 0; d < n_dg; d++)
        for (size_t gin = 0; gin < n_g; gin++) {
          tensor::View<double> row = chi_delayed.slice(a, d, gin);
          row /= row.sum();
        }
  }

  // Allow the ScattData object to combine itself
  for (size_t a = 0; a < total.shape(0); a++) {
    // Build vector of the scattering objects to incorporate
    vector<ScattData*> those_scatts(those_xs.size());
    for (size_t i = 0; i < those_xs.size(); i++) {
      those_scatts[i] = those_xs[i]->scatter[a].get();
    }

    // Now combine these guys
    scatter[a]->combine(those_scatts, scalars);
  }
}

//==============================================================================

bool XsData::equiv(const XsData& that)
{
  return (absorption.shape() == that.absorption.shape());
}

} // namespace openmc
