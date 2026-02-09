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

void XsData::from_hdf5(hid_t xsdata_grp, bool fissionable,
  AngleDistributionType scatter_format,
  AngleDistributionType final_scatter_format, int order_data, bool is_isotropic,
  int n_pol, int n_azi)
{
  // Reconstruct the dimension information so it doesn't need to be passed
  size_t n_ang = n_pol * n_azi;
  size_t energy_groups = total.shape()[1];

  // Set the fissionable-specific data
  if (fissionable) {
    fission_from_hdf5(xsdata_grp, n_ang, is_isotropic);
  }
  // Get the non-fission-specific data
  read_nd_vector(xsdata_grp, "decay-rate", decay_rate);
  read_nd_vector(xsdata_grp, "absorption", absorption, true);
  read_nd_vector(xsdata_grp, "inverse-velocity", inverse_velocity);

  // Get scattering data
  scatter_from_hdf5(
    xsdata_grp, n_ang, scatter_format, final_scatter_format, order_data);

  // Check absorption to ensure it is not 0 since it is often the
  // denominator in tally methods
  for (size_t i = 0; i < absorption.size(); ++i)
    if (absorption.data()[i] == 0.)
      absorption.data()[i] = 1.e-10;

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
  for (size_t i = 0; i < total.size(); ++i)
    if (total.data()[i] == 0.)
      total.data()[i] = 1.e-10;
}

//==============================================================================

void XsData::fission_vector_beta_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
{
  // Data is provided as nu-fission and chi with a beta for delayed info

  // Get chi
  xt::xtensor<double, 2> temp_chi({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  for (size_t a = 0; a < n_ang; ++a) {
    double row_sum = 0.0;
    for (size_t g = 0; g < n_g_; ++g)
      row_sum += temp_chi(a, g);
    for (size_t g = 0; g < n_g_; ++g)
      temp_chi(a, g) /= row_sum;
  }

  // Now every incoming group in prompt_chi and delayed_chi is the normalized
  // chi we just made (broadcast 2D -> 3D and 2D -> 4D)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin)
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) = temp_chi(a, gout);
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t d = 0; d < n_dg_; ++d)
      for (size_t gin = 0; gin < n_g_; ++gin)
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_delayed(a, d, gin, gout) = temp_chi(a, gout);

  // Get nu-fission
  xt::xtensor<double, 2> temp_nufiss({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "nu-fission", temp_nufiss, true);

  // Get beta (strategy will depend upon the number of dimensions in beta)
  hid_t beta_dset = open_dataset(xsdata_grp, "beta");
  int beta_ndims = dataset_ndims(beta_dset);
  close_dataset(beta_dset);
  int ndim_target = 1;
  if (!is_isotropic)
    ndim_target += 2;
  if (beta_ndims == ndim_target) {
    xt::xtensor<double, 2> temp_beta({n_ang, n_dg_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // Set prompt_nu_fission = (1. - beta_total)*nu_fission
    // beta_total = sum over delayed groups (axis 1) of temp_beta
    for (size_t a = 0; a < n_ang; ++a) {
      double beta_total = 0.0;
      for (size_t d = 0; d < n_dg_; ++d)
        beta_total += temp_beta(a, d);
      for (size_t g = 0; g < n_g_; ++g)
        prompt_nu_fission(a, g) = temp_nufiss(a, g) * (1.0 - beta_total);
    }

    // Set delayed_nu_fission as beta * nu_fission
    // delayed_nu_fission(a, d, g) = temp_beta(a, d) * temp_nufiss(a, g)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t g = 0; g < n_g_; ++g)
          delayed_nu_fission(a, d, g) = temp_beta(a, d) * temp_nufiss(a, g);
  } else if (beta_ndims == ndim_target + 1) {
    xt::xtensor<double, 3> temp_beta({n_ang, n_dg_, n_g_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // Set prompt_nu_fission = (1. - beta_total)*nu_fission
    // beta_total = sum over delayed groups (axis 1) of 3D temp_beta
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t g = 0; g < n_g_; ++g) {
        double beta_total = 0.0;
        for (size_t d = 0; d < n_dg_; ++d)
          beta_total += temp_beta(a, d, g);
        prompt_nu_fission(a, g) = temp_nufiss(a, g) * (1.0 - beta_total);
      }

    // Set delayed_nu_fission as beta * nu_fission
    // delayed_nu_fission(a, d, g) = temp_beta(a, d, g) * temp_nufiss(a, g)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t g = 0; g < n_g_; ++g)
          delayed_nu_fission(a, d, g) = temp_beta(a, d, g) * temp_nufiss(a, g);
  }
}

void XsData::fission_vector_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get chi-prompt
  xt::xtensor<double, 2> temp_chi_p({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi-prompt", temp_chi_p, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  for (size_t a = 0; a < n_ang; ++a) {
    double row_sum = 0.0;
    for (size_t g = 0; g < n_g_; ++g)
      row_sum += temp_chi_p(a, g);
    for (size_t g = 0; g < n_g_; ++g)
      temp_chi_p(a, g) /= row_sum;
  }

  // Get chi-delayed
  xt::xtensor<double, 3> temp_chi_d({n_ang, n_dg_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi-delayed", temp_chi_d, true);

  // Normalize chi by summing over the outgoing groups for each delayed group
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t d = 0; d < n_dg_; ++d) {
      double grp_sum = 0.0;
      for (size_t g = 0; g < n_g_; ++g)
        grp_sum += temp_chi_d(a, d, g);
      for (size_t g = 0; g < n_g_; ++g)
        temp_chi_d(a, d, g) /= grp_sum;
    }

  // Now assign the prompt and delayed chis by replicating for each incoming
  // group (broadcast 2D -> 3D and 3D -> 4D)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin)
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) = temp_chi_p(a, gout);
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t d = 0; d < n_dg_; ++d)
      for (size_t gin = 0; gin < n_g_; ++gin)
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_delayed(a, d, gin, gout) = temp_chi_d(a, d, gout);

  // Get prompt and delayed nu-fission directly
  read_nd_vector(xsdata_grp, "prompt-nu-fission", prompt_nu_fission, true);
  read_nd_vector(xsdata_grp, "delayed-nu-fission", delayed_nu_fission, true);
}

void XsData::fission_vector_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get chi
  xt::xtensor<double, 2> temp_chi({n_ang, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "chi", temp_chi, true);

  // Normalize chi by summing over the outgoing groups for each incoming angle
  for (size_t a = 0; a < n_ang; ++a) {
    double row_sum = 0.0;
    for (size_t g = 0; g < n_g_; ++g)
      row_sum += temp_chi(a, g);
    for (size_t g = 0; g < n_g_; ++g)
      temp_chi(a, g) /= row_sum;
  }

  // Now every incoming group in self.chi is the normalized chi we just made
  // (broadcast 2D -> 3D)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin)
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) = temp_chi(a, gout);

  // Get nu-fission directly
  read_nd_vector(xsdata_grp, "nu-fission", prompt_nu_fission, true);
}

//==============================================================================

void XsData::fission_matrix_beta_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
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
  if (!is_isotropic)
    ndim_target += 2;
  if (beta_ndims == ndim_target) {
    xt::xtensor<double, 2> temp_beta({n_ang, n_dg_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // temp_beta_sum(a) = sum over delayed groups of temp_beta(a, d)
    auto temp_beta_sum = xt::sum(temp_beta, {1});
    // matrix_sum(a, gin) = sum over outgoing groups of temp_matrix(a, gin, gout)
    auto matrix_sum = xt::sum(temp_matrix, {2});

    // prompt_nu_fission(a, g) = matrix_sum(a, g) * (1 - beta_sum(a))
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t g = 0; g < n_g_; ++g)
        prompt_nu_fission(a, g) = matrix_sum(a, g) * (1.0 - temp_beta_sum(a));

    // chi_prompt(a, gin, gout) = (1 - beta_sum(a)) * matrix(a, gin, gout)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t gin = 0; gin < n_g_; ++gin)
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_prompt(a, gin, gout) =
            (1.0 - temp_beta_sum(a)) * temp_matrix(a, gin, gout);

    // delayed_nu_fission(a, d, g) = beta(a, d) * matrix_sum(a, g)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t g = 0; g < n_g_; ++g)
          delayed_nu_fission(a, d, g) = temp_beta(a, d) * matrix_sum(a, g);

    // chi_delayed(a, d, gin, gout) = beta(a, d) * matrix(a, gin, gout)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t gin = 0; gin < n_g_; ++gin)
          for (size_t gout = 0; gout < n_g_; ++gout)
            chi_delayed(a, d, gin, gout) =
              temp_beta(a, d) * temp_matrix(a, gin, gout);

  } else if (beta_ndims == ndim_target + 1) {
    xt::xtensor<double, 3> temp_beta({n_ang, n_dg_, n_g_}, 0.);
    read_nd_vector(xsdata_grp, "beta", temp_beta, true);

    // temp_beta_sum(a, g) = sum over delayed groups of temp_beta(a, d, g)
    auto temp_beta_sum = xt::sum(temp_beta, {1});
    auto matrix_sum = xt::sum(temp_matrix, {2});

    // prompt_nu_fission(a, g) = matrix_sum(a, g) * (1 - beta_sum(a, g))
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t g = 0; g < n_g_; ++g)
        prompt_nu_fission(a, g) = matrix_sum(a, g) * (1.0 - temp_beta_sum(a, g));

    // chi_prompt(a, gin, gout) = (1 - beta_sum(a, gin)) * matrix(a, gin, gout)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t gin = 0; gin < n_g_; ++gin)
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_prompt(a, gin, gout) =
            (1.0 - temp_beta_sum(a, gin)) * temp_matrix(a, gin, gout);

    // delayed_nu_fission(a, d, g) = beta(a, d, g) * matrix_sum(a, g)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t g = 0; g < n_g_; ++g)
          delayed_nu_fission(a, d, g) = temp_beta(a, d, g) * matrix_sum(a, g);

    // chi_delayed(a, d, gin, gout) = beta(a, d, gin) * matrix(a, gin, gout)
    for (size_t a = 0; a < n_ang; ++a)
      for (size_t d = 0; d < n_dg_; ++d)
        for (size_t gin = 0; gin < n_g_; ++gin)
          for (size_t gout = 0; gout < n_g_; ++gout)
            chi_delayed(a, d, gin, gout) =
              temp_beta(a, d, gin) * temp_matrix(a, gin, gout);
  }

  // Normalize both chis: chi_prompt(a, gin, gout) /= sum_over_gout
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin) {
      double s = 0.0;
      for (size_t gout = 0; gout < n_g_; ++gout)
        s += chi_prompt(a, gin, gout);
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) /= s;
    }

  // chi_delayed(a, d, gin, gout) /= sum_over_gout
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t d = 0; d < n_dg_; ++d)
      for (size_t gin = 0; gin < n_g_; ++gin) {
        double s = 0.0;
        for (size_t gout = 0; gout < n_g_; ++gout)
          s += chi_delayed(a, d, gin, gout);
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_delayed(a, d, gin, gout) /= s;
      }
}

void XsData::fission_matrix_no_beta_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // Data is provided separately as prompt + delayed nu-fission and chi

  // Get the prompt nu-fission matrix
  xt::xtensor<double, 3> temp_matrix_p({n_ang, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "prompt-nu-fission", temp_matrix_p, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = xt::sum(temp_matrix_p, {2});

  // chi_prompt = matrix / prompt_nu_fission (broadcast over gout)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin)
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) =
          temp_matrix_p(a, gin, gout) / prompt_nu_fission(a, gin);

  // Get the delayed nu-fission matrix
  xt::xtensor<double, 4> temp_matrix_d({n_ang, n_dg_, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_matrix_d, true);

  // delayed_nu_fission is the sum over outgoing groups
  delayed_nu_fission = xt::sum(temp_matrix_d, {3});

  // chi_delayed = matrix / delayed_nu_fission (broadcast over gout)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t d = 0; d < n_dg_; ++d)
      for (size_t gin = 0; gin < n_g_; ++gin)
        for (size_t gout = 0; gout < n_g_; ++gout)
          chi_delayed(a, d, gin, gout) =
            temp_matrix_d(a, d, gin, gout) / delayed_nu_fission(a, d, gin);
}

void XsData::fission_matrix_no_delayed_from_hdf5(hid_t xsdata_grp, size_t n_ang)
{
  // No beta is provided and there is no prompt/delay distinction.
  // Therefore, the code only considers the data as prompt.

  // Get nu-fission matrix
  xt::xtensor<double, 3> temp_matrix({n_ang, n_g_, n_g_}, 0.);
  read_nd_vector(xsdata_grp, "nu-fission", temp_matrix, true);

  // prompt_nu_fission is the sum over outgoing groups
  prompt_nu_fission = xt::sum(temp_matrix, {2});

  // chi_prompt = matrix / prompt_nu_fission (broadcast over gout)
  for (size_t a = 0; a < n_ang; ++a)
    for (size_t gin = 0; gin < n_g_; ++gin)
      for (size_t gout = 0; gout < n_g_; ++gout)
        chi_prompt(a, gin, gout) =
          temp_matrix(a, gin, gout) / prompt_nu_fission(a, gin);
}

//==============================================================================

void XsData::fission_from_hdf5(
  hid_t xsdata_grp, size_t n_ang, bool is_isotropic)
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
    // nu_fission = prompt_nu_fission + sum over delayed groups
    nu_fission = prompt_nu_fission;
    for (size_t a = 0; a < nu_fission.shape(0); ++a)
      for (size_t g = 0; g < nu_fission.shape(1); ++g)
        for (size_t d = 0; d < delayed_nu_fission.shape(1); ++d)
          nu_fission(a, g) += delayed_nu_fission(a, d, g);
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
  xt::xtensor<int, 2> gmin({n_ang, n_g_}, 0.);
  read_nd_vector(scatt_grp, "g_min", gmin, true);
  xt::xtensor<int, 2> gmax({n_ang, n_g_}, 0.);
  read_nd_vector(scatt_grp, "g_max", gmax, true);

  // Make gmin and gmax start from 0 vice 1 as they do in the library
  gmin -= 1;
  gmax -= 1;

  // Now use this info to find the length of a vector to hold the flattened
  // data.
  size_t length = order_data * static_cast<size_t>(xt::sum(gmax - gmin + 1));

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
      xt::xtensor<int, 1> in_gmin = xt::view(gmin, a, xt::all());
      xt::xtensor<int, 1> in_gmax = xt::view(gmax, a, xt::all());
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
    if (that->prompt_nu_fission.shape()[0] > 0) {
      nu_fission += scalar * that->nu_fission;
      prompt_nu_fission += scalar * that->prompt_nu_fission;
      kappa_fission += scalar * that->kappa_fission;
      fission += scalar * that->fission;
      delayed_nu_fission += scalar * that->delayed_nu_fission;
      // chi_prompt += scalar * sum_pnf(a) * that->chi_prompt(a, gin, gout)
      // where sum_pnf(a) = sum over gin of prompt_nu_fission(a, gin)
      {
        auto sum_pnf = xt::sum(that->prompt_nu_fission, {1});
        for (size_t a = 0; a < chi_prompt.shape(0); ++a)
          for (size_t gin = 0; gin < chi_prompt.shape(1); ++gin)
            for (size_t gout = 0; gout < chi_prompt.shape(2); ++gout)
              chi_prompt(a, gin, gout) +=
                scalar * sum_pnf(a) * that->chi_prompt(a, gin, gout);
      }
      // chi_delayed += scalar * sum_dnf(a, d) * that->chi_delayed(a, d, gin, gout)
      // where sum_dnf(a, d) = sum over gin of delayed_nu_fission(a, d, gin)
      {
        auto sum_dnf = xt::sum(that->delayed_nu_fission, {2});
        for (size_t a = 0; a < chi_delayed.shape(0); ++a)
          for (size_t d = 0; d < chi_delayed.shape(1); ++d)
            for (size_t gin = 0; gin < chi_delayed.shape(2); ++gin)
              for (size_t gout = 0; gout < chi_delayed.shape(3); ++gout)
                chi_delayed(a, d, gin, gout) +=
                  scalar * sum_dnf(a, d) * that->chi_delayed(a, d, gin, gout);
      }
    }
    decay_rate += scalar * that->decay_rate;
  }

  // Ensure the chi_prompt and chi_delayed are normalized to 1
  for (size_t a = 0; a < chi_prompt.shape(0); ++a)
    for (size_t gin = 0; gin < chi_prompt.shape(1); ++gin) {
      double s = 0.0;
      for (size_t gout = 0; gout < chi_prompt.shape(2); ++gout)
        s += chi_prompt(a, gin, gout);
      for (size_t gout = 0; gout < chi_prompt.shape(2); ++gout)
        chi_prompt(a, gin, gout) /= s;
    }
  for (size_t a = 0; a < chi_delayed.shape(0); ++a)
    for (size_t d = 0; d < chi_delayed.shape(1); ++d)
      for (size_t gin = 0; gin < chi_delayed.shape(2); ++gin) {
        double s = 0.0;
        for (size_t gout = 0; gout < chi_delayed.shape(3); ++gout)
          s += chi_delayed(a, d, gin, gout);
        for (size_t gout = 0; gout < chi_delayed.shape(3); ++gout)
          chi_delayed(a, d, gin, gout) /= s;
      }

  // Allow the ScattData object to combine itself
  for (size_t a = 0; a < total.shape()[0]; a++) {
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
