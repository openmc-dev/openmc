#include "openmc/xsdata.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"


namespace openmc {

//==============================================================================
// XsData class methods
//==============================================================================

XsData::XsData(int energy_groups, int num_delayed_groups, bool fissionable,
     int scatter_format, int n_pol, int n_azi)
{
  size_t n_ang = n_pol * n_azi;

  // check to make sure scatter format is OK before we allocate
  if (scatter_format != ANGLE_HISTOGRAM && scatter_format != ANGLE_TABULAR &&
      scatter_format != ANGLE_LEGENDRE) {
    fatal_error("Invalid scatter_format!");
  }
  // allocate all [temperature][phi][theta][in group] quantities
  std::vector<size_t> shape = {n_ang, energy_groups};
  total = xt::zeros<double>(shape);
  absorption = xt::zeros<double>(shape);
  inverse_velocity = xt::zeros<double>(shape);
  if (fissionable) {
    fission = xt::zeros<double>(shape);
    nu_fission = xt::zeros<double>(shape);
    prompt_nu_fission = xt::zeros<double>(shape);
    kappa_fission = xt::zeros<double>(shape);
  }

  // allocate decay_rate; [temperature][phi][theta][delayed group]
  shape[1] = num_delayed_groups;
  decay_rate = xt::zeros<double>(shape);

  if (fissionable) {
    shape = {n_ang, energy_groups, num_delayed_groups};
    // allocate delayed_nu_fission; [temperature][phi][theta][in group][delay group]
    delayed_nu_fission = xt::zeros<double>(shape);

    // chi_prompt; [temperature][phi][theta][in group][out group]
    shape = {n_ang, energy_groups, energy_groups};
    chi_prompt = xt::zeros<double>(shape);

    // chi_delayed; [temperature][phi][theta][in group][out group][delay group]
    shape = {n_ang, energy_groups, energy_groups, num_delayed_groups};
    chi_delayed = xt::zeros<double>(shape);
  }


  for (int a = 0; a < n_ang; a++) {
    if (scatter_format == ANGLE_HISTOGRAM) {
      // scatter[a] = std::make_unique(ScattDataHistogram);
      scatter.emplace_back(new ScattDataHistogram);
    } else if (scatter_format == ANGLE_TABULAR) {
      // scatter[a] = std::make_unique(ScattDataTabular);
      scatter.emplace_back(new ScattDataTabular);
    } else if (scatter_format == ANGLE_LEGENDRE) {
      // scatter[a] = std::make_unique(ScattDataLegendre);
      scatter.emplace_back(new ScattDataLegendre);
    }
  }
}

//==============================================================================

void
XsData::from_hdf5(hid_t xsdata_grp, bool fissionable, int scatter_format,
     int final_scatter_format, int order_data, int max_order,
     int legendre_to_tabular_points, bool is_isotropic, int n_pol, int n_azi)
{
  // Reconstruct the dimension information so it doesn't need to be passed
  int n_ang = n_pol * n_azi;
  int energy_groups = total.shape()[1];
  int delayed_groups = decay_rate.shape()[1];

  // Set the fissionable-specific data
  if (fissionable) {
    fission_from_hdf5(xsdata_grp, n_pol, n_azi, energy_groups, delayed_groups,
                      is_isotropic);
  }
  // Get the non-fission-specific data
  read_nd_vector(xsdata_grp, "decay_rate", decay_rate);
  read_nd_vector(xsdata_grp, "absorption", absorption, true);
  read_nd_vector(xsdata_grp, "inverse-velocity", inverse_velocity);

  // Get scattering data
  scatter_from_hdf5(xsdata_grp, n_pol, n_azi, energy_groups, scatter_format,
       final_scatter_format, order_data, max_order, legendre_to_tabular_points);

  // Check absorption to ensure it is not 0 since it is often the
  // denominator in tally methods
  for (int a = 0; a < n_ang; a++) {
    for (int gin = 0; gin < energy_groups; gin++) {
      if (absorption(a, gin) == 0.) absorption(a, gin) = 1.e-10;
    }
  }

  // Get or calculate the total x/s
  if (object_exists(xsdata_grp, "total")) {
    read_nd_vector(xsdata_grp, "total", total);
  } else {
    for (int a = 0; a < n_ang; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        total(a, gin) = absorption(a, gin) + scatter[a]->scattxs[gin];
      }
    }
  }

  // Fix if total is 0, since it is in the denominator when tallying
  for (int a = 0; a < n_ang; a++) {
    for (int gin = 0; gin < energy_groups; gin++) {
      if (total(a, gin) == 0.) total(a, gin) = 1.e-10;
    }
  }
}

//==============================================================================

void
XsData::fission_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
     int energy_groups, int delayed_groups, bool is_isotropic)
{
  size_t n_ang = n_pol * n_azi;
  // Get the fission and kappa_fission data xs; these are optional
  read_nd_vector(xsdata_grp, "fission", fission);
  read_nd_vector(xsdata_grp, "kappa-fission", kappa_fission);

  // Set/get beta
  xt::xtensor<double, 3> temp_beta({n_ang, energy_groups, delayed_groups}, 0.);
  if (object_exists(xsdata_grp, "beta")) {
    hid_t xsdata = open_dataset(xsdata_grp, "beta");
    int ndims = dataset_ndims(xsdata);

    // raise ndims to make the isotropic ndims the same as angular
    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // Beta is input as [delayed group]
      std::vector<double> temp_arr({n_pol * n_azi * delayed_groups});
      read_nd_vector(xsdata_grp, "beta", temp_arr);

      // Broadcast to all incoming groups
      int temp_idx = 0;
      for (int a = 0; a < n_ang; a++) {
        for (int dg = 0; dg < delayed_groups; dg++) {
          // Set the first group index and copy the rest
          temp_beta(a, 0, dg) = temp_arr[temp_idx++];
          for (int gin = 1; gin < energy_groups; gin++) {
            temp_beta(a, gin, dg) = temp_beta(a, 0, dg);
          }
        }
      }
    } else if (ndims == 4) {
      // Beta is input as [in group][delayed group]
      read_nd_vector(xsdata_grp, "beta", temp_beta);
    } else {
      fatal_error("beta must be provided as a 3D or 4D array!");
    }
  }

  // If chi is provided, set chi-prompt and chi-delayed
  if (object_exists(xsdata_grp, "chi")) {
    xt::xtensor<double, 2> temp_arr ({n_ang, energy_groups});
    read_nd_vector(xsdata_grp, "chi", temp_arr);

    for (int a = 0; a < n_ang; a++) {
      // First set the first group
      for (int gout = 0; gout < energy_groups; gout++) {
        chi_prompt(a, 0, gout) = temp_arr(a, gout);
      }

      // Now normalize this data
      double chi_sum = 0.;
      for (int gout = 0; gout < energy_groups; gout++) {
        chi_sum += chi_prompt(a, 0, gout);
      }

      if (chi_sum <= 0.) {
        fatal_error("Encountered chi for a group that is <= 0!");
      }
      for (int gout = 0; gout < energy_groups; gout++) {
        chi_prompt(a, 0, gout) /= chi_sum;
      }

      // And extend to the remaining incoming groups
      for (int gin = 1; gin < energy_groups; gin++) {
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_prompt(a, gin, gout) = chi_prompt(a, 0, gout);
        }
      }

      // Finally set chi-delayed equal to chi-prompt
      // Set chi-delayed to chi-prompt
      for(int gin = 0; gin < energy_groups; gin++) {
        for (int gout = 0; gout < energy_groups; gout++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            chi_delayed(a, gin, gout, dg) = chi_prompt(a, gin, gout);
          }
        }
      }
    }
  }

  // If nu-fission is provided, set prompt- and delayed-nu-fission;
  // if nu-fission is a matrix, set chi-prompt and chi-delayed.
  if (object_exists(xsdata_grp, "nu-fission")) {
    hid_t xsdata = open_dataset(xsdata_grp, "nu-fission");
    int ndims = dataset_ndims(xsdata);
    // raise ndims to make the isotropic ndims the same as angular
    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // nu-fission is a 3-d array
      read_nd_vector(xsdata_grp, "nu-fission", prompt_nu_fission);

      // set delayed-nu-fission and correct prompt-nu-fission with beta
      for (int a = 0; a < n_ang; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            delayed_nu_fission(a, gin, dg) =
                 temp_beta(a, gin, dg) * prompt_nu_fission(a, gin);
          }

          // Correct the prompt-nu-fission using the delayed neutron fraction
          if (delayed_groups > 0) {
            double beta_sum = 0.;
            for (int gin = 0; gin < energy_groups; gin++) {
              beta_sum += temp_beta(a, gin);
            }

            prompt_nu_fission(a, gin) *= (1. - beta_sum);
          }
        }
      }

    } else if (ndims == 4) {
      // nu-fission is a matrix
      read_nd_vector(xsdata_grp, "nu-fission", chi_prompt);

      // Normalize the chi info so the CDF is 1.
      for (int a = 0; a < n_ang; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          double chi_sum = 0.;
          for (int gout = 0; gout < energy_groups; gout++) {
            chi_sum += chi_prompt(a, gin, gout);
          }

          // Set the vector nu-fission from the matrix nu-fission
          prompt_nu_fission(a, gin) = chi_sum;

          if (chi_sum >= 0.) {
            for (int gout = 0; gout < energy_groups; gout++) {
              chi_prompt(a, gin, gout) /= chi_sum;
            }
          } else {
            fatal_error("Encountered chi for a group that is <= 0!");
          }
        }

        // set chi-delayed to chi-prompt
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int gout = 0; gout < energy_groups; gout++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              chi_delayed(a, gin, gout, dg) = chi_prompt(a, gin, gout);
            }
          }
        }

        // Set the delayed-nu-fission and correct prompt-nu-fission with beta
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            delayed_nu_fission(a, gin, dg) = temp_beta(a, gin, dg) *
                 prompt_nu_fission(a, gin);
          }

          // Correct prompt-nu-fission using the delayed neutron fraction
          if (delayed_groups > 0) {
            double beta_sum = 0.;
            for (int dg = 0; dg < delayed_groups; dg++) {
              beta_sum += temp_beta(a, gin, dg);
            }
            prompt_nu_fission(a, gin) *= (1. - beta_sum);
          }
        }
      }
    } else {
      fatal_error("nu-fission must be provided as a 3D or 4D array!");
    }

    close_dataset(xsdata);
  }

  // If chi-prompt is provided, set chi-prompt
  if (object_exists(xsdata_grp, "chi-prompt")) {
    xt::xtensor<double, 2> temp_arr({n_ang, energy_groups});
    read_nd_vector(xsdata_grp, "chi-prompt", temp_arr);

    for (int a = 0; a < n_ang; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_prompt(a, gin, gout) = temp_arr(a, gout);
        }

        // Normalize chi so its CDF goes to 1
        double chi_sum = 0.;
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_sum += chi_prompt(a, gin, gout);
        }

        if (chi_sum >= 0.) {
          for (int gout = 0; gout < energy_groups; gout++) {
            chi_prompt(a, gin, gout) /= chi_sum;
          }
        } else {
          fatal_error("Encountered chi-prompt for a group that is <= 0.!");
        }
      }
    }
  }

  // If chi-delayed is provided, set chi-delayed
  if (object_exists(xsdata_grp, "chi-delayed")) {
    hid_t xsdata = open_dataset(xsdata_grp, "chi-delayed");
    int ndims = dataset_ndims(xsdata);
    // raise ndims to make the isotropic ndims the same as angular
    if (is_isotropic) ndims += 2;
    close_dataset(xsdata);

    if (ndims == 3) {
      // chi-delayed is a [in group] vector
      xt::xtensor<double, 2> temp_arr({n_ang, energy_groups});
      read_nd_vector(xsdata_grp, "chi-delayed", temp_arr);

      for (int a = 0; a < n_ang; a++) {
        // normalize the chi CDF to 1
        double chi_sum = 0.;
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_sum += temp_arr(a, gout);
        }

        if (chi_sum <= 0.) {
          fatal_error("Encountered chi-delayed for a group that is <= 0!");
        }

        // set chi-delayed
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int gout = 0; gout < energy_groups; gout++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              chi_delayed(a, gin, gout, dg) = temp_arr(a, gout) / chi_sum;
            }
          }
        }
      }
    } else if (ndims == 4) {
      // chi_delayed is a matrix
      read_nd_vector(xsdata_grp, "chi-delayed", chi_delayed);

      // Normalize the chi info so the CDF is 1.
      for (int a = 0; a < n_ang; a++) {
        for (int dg = 0; dg < delayed_groups; dg++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            double chi_sum = 0.;
            for (int gout = 0; gout < energy_groups; gout++) {
              chi_sum += chi_delayed(a, gin, gout, dg);
            }

            if (chi_sum > 0.) {
              for (int gout = 0; gout < energy_groups; gout++) {
                chi_delayed(a, gin, gout, dg) /= chi_sum;
              }
            } else {
              fatal_error("Encountered chi-delayed for a group that is <= 0!");
            }
          }
        }
      }
    } else {
      fatal_error("chi-delayed must be provided as a 3D or 4D array!");
    }
  }

  // Get prompt-nu-fission, if present
  if (object_exists(xsdata_grp, "prompt-nu-fission")) {
    hid_t xsdata = open_dataset(xsdata_grp, "prompt-nu-fission");
    int ndims = dataset_ndims(xsdata);
    // raise ndims to make the isotropic ndims the same as angular
    if (is_isotropic) ndims += 2;
    close_dataset(xsdata);

    if (ndims == 3) {
      // prompt-nu-fission is a [in group] vector
      read_nd_vector(xsdata_grp, "prompt-nu-fission", prompt_nu_fission);
    } else if (ndims == 4) {
      // prompt nu fission is a matrix,
      // so set prompt_nu_fiss & chi_prompt
      xt::xtensor<double, 3> temp_arr({n_ang, energy_groups, energy_groups});
      read_nd_vector(xsdata_grp, "prompt-nu-fission", temp_arr);

      // The prompt_nu_fission vector from the matrix form
      for (int a = 0; a < n_ang; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          double prompt_sum = 0.;
          for (int gout = 0; gout < energy_groups; gout++) {
            prompt_sum += temp_arr(a, gin, gout);
          }

          prompt_nu_fission(a, gin) = prompt_sum;
        }

        // The chi_prompt data is just the normalized fission matrix
        for (int gin= 0; gin < energy_groups; gin++) {
          if (prompt_nu_fission(a, gin) > 0.) {
            for (int gout = 0; gout < energy_groups; gout++) {
              chi_prompt(a, gin, gout) =
                   temp_arr(a, gin, gout) / prompt_nu_fission(a, gin);
            }
          } else {
            fatal_error("Encountered chi-prompt for a group that is <= 0!");
          }
        }
      }

    } else {
      fatal_error("prompt-nu-fission must be provided as a 3D or 4D array!");
    }
  }

  // Get delayed-nu-fission, if present
  if (object_exists(xsdata_grp, "delayed-nu-fission")) {
    hid_t xsdata = open_dataset(xsdata_grp, "delayed-nu-fission");
    int ndims = dataset_ndims(xsdata);
    close_dataset(xsdata);
    // raise ndims to make the isotropic ndims the same as angular
    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // delayed-nu-fission is an [in group] vector
      if (temp_beta(0, 0, 0) == 0.) {
        fatal_error("cannot set delayed-nu-fission with a 1D array if "
                    "beta is not provided");
      }
      xt::xtensor<double, 2> temp_arr({n_ang, energy_groups});
      read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_arr);

      for (int a = 0; a < n_ang; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            // Set delayed-nu-fission using beta
            delayed_nu_fission(a, gin, dg) =
                 temp_beta(a, gin, dg) * temp_arr(a, gin);
          }
        }
      }

    } else if (ndims == 4) {
      read_nd_vector(xsdata_grp, "delayed-nu-fission",
                    delayed_nu_fission);

    } else if (ndims == 5) {
      // This will contain delayed-nu-fission and chi-delayed data
      xt::xtensor<double, 4> temp_arr({n_ang, energy_groups, energy_groups,
                                       delayed_groups});
      read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_arr);

      // Set the 3D delayed-nu-fission matrix and 4D chi-delayed matrix
      // from the 4D delayed-nu-fission matrix
      for (int a = 0; a < n_ang; a++) {
        for (int dg = 0; dg < delayed_groups; dg++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            double gout_sum = 0.;
            for (int gout = 0; gout < energy_groups; gout++) {
              gout_sum += temp_arr(a, gin, gout, dg);
              chi_delayed(a, gin, gout, dg) = temp_arr(a, gin, gout, dg);
            }
            delayed_nu_fission(a, gin, dg) = gout_sum;
            // Normalize chi-delayed
            if (gout_sum > 0.) {
              for (int gout = 0; gout < energy_groups; gout++) {
                chi_delayed(a, gin, gout, dg) /= gout_sum;
              }
            } else {
              fatal_error("Encountered chi-delayed for a group that is <= 0!");
            }
          }
        }
      }

    } else {
      fatal_error("prompt-nu-fission must be provided as a 3D, 4D, or 5D "
                  "array!");
    }
  }

  // Combine prompt_nu_fission and delayed_nu_fission into nu_fission
  for (int a = 0; a < n_ang; a++) {
    for (int gin = 0; gin < energy_groups; gin++) {
      nu_fission(a, gin) = prompt_nu_fission(a, gin);
      for (int dg = 0; dg < delayed_groups; dg++) {
        nu_fission(a, gin) += delayed_nu_fission(a, gin, dg);
      }
    }
  }
}

//==============================================================================

void
XsData::scatter_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
     int energy_groups, int scatter_format, int final_scatter_format,
     int order_data, int max_order, int legendre_to_tabular_points)
{
  size_t n_ang = n_pol * n_azi;
  if (!object_exists(xsdata_grp, "scatter_data")) {
    fatal_error("Must provide scatter_data group!");
  }
  hid_t scatt_grp = open_group(xsdata_grp, "scatter_data");

  // Get the outgoing group boundary indices
  xt::xtensor<int, 2> gmin({n_ang, energy_groups});
  read_nd_vector(scatt_grp, "g_min", gmin, true);
  xt::xtensor<int, 2> gmax({n_ang, energy_groups});
  read_nd_vector(scatt_grp, "g_max", gmax, true);

  // Make gmin and gmax start from 0 vice 1 as they do in the library
  gmin -= 1;
  gmax -= 1;

  // Now use this info to find the length of a vector to hold the flattened
  // data.
  int length = 0;
  for (int a = 0; a < n_ang; a++) {
    for (int gin = 0; gin < energy_groups; gin++) {
      length += order_data * (gmax(a, gin) - gmin(a, gin) + 1);
    }
  }

  double_4dvec input_scatt(n_ang, double_3dvec(energy_groups));
  //temp_arr scope
  {
    std::vector<double> temp_arr(length);
    read_nd_vector(scatt_grp, "scatter_matrix", temp_arr, true);

    // Compare the number of orders given with the max order of the problem;
    // strip off the superfluous orders if needed
    int order_dim;
    if (scatter_format == ANGLE_LEGENDRE) {
      order_dim = std::min(order_data - 1, max_order) + 1;
    } else {
      order_dim = order_data;
    }

    // convert the flattened temp_arr to a jagged array for passing to
    // scatt data
    int temp_idx = 0;
    for (int a = 0; a < n_ang; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        input_scatt[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
        for (int i_gout = 0; i_gout < input_scatt[a][gin].size(); i_gout++) {
          input_scatt[a][gin][i_gout].resize(order_dim);
          for (int l = 0; l < order_dim; l++) {
            input_scatt[a][gin][i_gout][l] = temp_arr[temp_idx++];
          }
          // Adjust index for the orders we didnt take
          temp_idx += (order_data - order_dim);
        }
      }
    }
  }

  // Get multiplication matrix
  double_3dvec temp_mult(n_ang, double_2dvec(energy_groups));
  if (object_exists(scatt_grp, "multiplicity_matrix")) {
    std::vector<double> temp_arr(length / order_data);
    read_nd_vector(scatt_grp, "multiplicity_matrix", temp_arr);

    // convert the flat temp_arr to a jagged array for passing to scatt data
    int temp_idx = 0;
    for (int a = 0; a < n_ang; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        temp_mult[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
        for (int i_gout = 0; i_gout < temp_mult[a][gin].size(); i_gout++) {
          temp_mult[a][gin][i_gout] = temp_arr[temp_idx++];
        }
      }
    }
    temp_arr.clear();
  } else {
    // Use a default: multiplicities are 1.0.
    for (int a = 0; a < n_ang; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        temp_mult[a][gin].resize(gmax(a, gin) - gmin(a, gin) + 1);
        for (int i_gout = 0; i_gout < temp_mult[a][gin].size(); i_gout++) {
          temp_mult[a][gin][i_gout] = 1.;
        }
      }
    }
  }
  close_group(scatt_grp);

  // Finally, convert the Legendre data to tabular, if needed
  if (scatter_format == ANGLE_LEGENDRE &&
      final_scatter_format == ANGLE_TABULAR) {
    for (int a = 0; a < n_ang; a++) {
      ScattDataLegendre legendre_scatt;
      xt::xtensor<int, 1> in_gmin = xt::view(gmin, a, xt::all());
      xt::xtensor<int, 1> in_gmax = xt::view(gmax, a, xt::all());
      legendre_scatt.init(in_gmin, in_gmax,
                          temp_mult[a], input_scatt[a]);

      // Now create a tabular version of legendre_scatt
      convert_legendre_to_tabular(legendre_scatt,
           *static_cast<ScattDataTabular*>(scatter[a].get()),
           legendre_to_tabular_points);

      scatter_format = final_scatter_format;
    }
  } else {
    // We are sticking with the current representation
    // Initialize the ScattData object with this data
    for (int a = 0; a < n_ang; a++) {
      scatter[a]->init(xt::view(gmin, a, xt::all()),
                       xt::view(gmax, a, xt::all()),
                       temp_mult[a], input_scatt[a]);
    }
  }
}

//==============================================================================

void
XsData::combine(const std::vector<XsData*>& those_xs,
                const double_1dvec& scalars)
{
  // Combine the non-scattering data
  for (int i = 0; i < those_xs.size(); i++) {
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
  for (int a = 0; a < total.shape()[0]; a++) {
    // Build vector of the scattering objects to incorporate
    std::vector<ScattData*> those_scatts(those_xs.size());
    for (int i = 0; i < those_xs.size(); i++) {
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
