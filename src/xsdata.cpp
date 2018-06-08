#include "xsdata.h"

namespace openmc {

//==============================================================================
// XsData class methods
//==============================================================================

XsData::XsData(int energy_groups, int num_delayed_groups, bool fissionable,
               int scatter_format, int n_pol, int n_azi)
{
  // check to make sure scatter format is OK before we allocate
  if (scatter_format != ANGLE_HISTOGRAM && scatter_format != ANGLE_TABULAR &&
      scatter_format != ANGLE_LEGENDRE) {
    fatal_error("Invalid scatter_format!");
  }
  // allocate all [temperature][phi][theta][in group] quantities
  total = double_3dvec(n_pol, double_2dvec(n_azi,
       double_1dvec(energy_groups, 0.)));
  absorption = double_3dvec(n_pol, double_2dvec(n_azi,
       double_1dvec(energy_groups, 0.)));
  inverse_velocity = double_3dvec(n_pol,
       double_2dvec(n_azi, double_1dvec(energy_groups, 0.)));
  if (fissionable) {
    fission = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups, 0.)));
    prompt_nu_fission = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups, 0.)));
    kappa_fission = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups, 0.)));
  }

  // allocate decay_rate; [temperature][phi][theta][delayed group]
  decay_rate = double_3dvec(n_pol, double_2dvec(n_azi,
       double_1dvec(num_delayed_groups, 0.)));

  if (fissionable) {
    // allocate delayed_nu_fission; [temperature][phi][theta][in group][out group]
    delayed_nu_fission = double_4dvec(n_pol, double_3dvec(n_azi,
         double_2dvec(energy_groups, double_1dvec(energy_groups, 0.))));

    // chi_prompt; [temperature][phi][theta][in group][delayed group]
    chi_prompt = double_4dvec(n_pol, double_3dvec(n_azi,
         double_2dvec(energy_groups, double_1dvec(energy_groups, 0.))));

    // chi_delayed; [temperature][phi][theta][in group][out group][delay group]
    chi_delayed = double_5dvec(n_pol, double_4dvec(n_azi,
         double_3dvec(energy_groups, double_2dvec(energy_groups,
         double_1dvec(num_delayed_groups, 0.)))));
  }

  scatter.resize(n_pol);
  for (int p = 0; p < n_pol; p++) {
    scatter[p].resize(n_azi);
    for (int a = 0; a < n_azi; a++) {
      if (scatter_format == ANGLE_HISTOGRAM) {
        scatter[p][a] = new ScattDataHistogram;
      } else if (scatter_format == ANGLE_TABULAR) {
        scatter[p][a] = new ScattDataTabular;
      } else if (scatter_format == ANGLE_LEGENDRE) {
        scatter[p][a] = new ScattDataLegendre;
      }
    }
  }
}


void XsData::from_hdf5(hid_t xsdata_grp, bool fissionable, int scatter_format,
                       int final_scatter_format, int order_data, int max_order,
                       int legendre_to_tabular_points, bool is_isotropic)
{
  // Reconstruct the dimension information so it doesn't need to be passed
  int n_pol = total.size();
  int n_azi = total[0].size();
  int energy_groups = total[0][0].size();
  int delayed_groups = decay_rate[0][0].size();

  // Set the fissionable-specific data
  if (fissionable) {
    _fissionable_from_hdf5(xsdata_grp, n_pol, n_azi, energy_groups,
                           delayed_groups, is_isotropic);
  }
  // Get the non-fission-specific data
  read_nd_vector(xsdata_grp, "decay_rate", decay_rate);
  read_nd_vector(xsdata_grp, "absorption", absorption, true);
  read_nd_vector(xsdata_grp, "inverse-velocity", inverse_velocity);

  // Get scattering data
  _scatter_from_hdf5(xsdata_grp, n_pol, n_azi, energy_groups, scatter_format,
                     final_scatter_format, order_data, max_order,
                     legendre_to_tabular_points);

  // Check absorption to ensure it is not 0 since it is often the
  // denominator in tally methods
  for (int p = 0; p < n_pol; p++) {
    for (int a = 0; a < n_azi; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        if (absorption[p][a][gin] == 0.) absorption[p][a][gin] = 1.e-10;
      }
    }
  }

  // Get or calculate the total x/s
  if (object_exists(xsdata_grp, "total")) {
    read_nd_vector(xsdata_grp, "total", total);
  } else {
    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          total[p][a][gin] = absorption[p][a][gin] +
               scatter[p][a]->scattxs[gin];
        }
      }
    }
  }

  // Check total to ensure it is not 0 since it is often the denominator in
  // tally methods
  for (int p = 0; p < n_pol; p++) {
    for (int a = 0; a < n_azi; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        if (total[p][a][gin] == 0.) total[p][a][gin] = 1.e-10;
      }
    }
  }
}


void XsData::_fissionable_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
     int energy_groups, int delayed_groups, bool is_isotropic)
{
  double_4dvec temp_beta =
       double_4dvec(n_pol, double_3dvec(n_azi,
       double_2dvec(energy_groups, double_1dvec(delayed_groups, 0.))));

  // Set/get beta
  if (object_exists(xsdata_grp, "beta")) {
    hid_t xsdata = open_dataset(xsdata_grp, "beta");
    int ndims = dataset_ndims(xsdata);

    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // Beta is input as [delayed group]
      double_1dvec temp_arr = double_1dvec(n_pol * n_azi * delayed_groups);
      read_nd_vector(xsdata_grp, "beta", temp_arr);

      // Broadcast to all incoming groups
      int temp_idx = 0;
      for (int p = 0; p < n_pol; p++) {
        for (int a = 0; a < n_azi; a++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            // Set the first group index and copy the rest
            temp_beta[p][a][0][dg] = temp_arr[temp_idx++];
            for (int gin = 1; gin < energy_groups; gin++) {
              temp_beta[p][a][gin] = temp_beta[p][a][0];
            }
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
    double_3dvec temp_arr = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups)));
    read_nd_vector(xsdata_grp, "chi", temp_arr);

    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        // First set the first group
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_prompt[p][a][0][gout] = temp_arr[p][a][gout];
        }

        // Now normalize this data
        double chi_sum = std::accumulate(chi_prompt[p][a][0].begin(),
                                         chi_prompt[p][a][0].end(),
                                         0.);
        if (chi_sum <= 0.) {
          fatal_error("Encountered chi for a group that is <= 0!");
        }
        for (int gout = 0; gout < energy_groups; gout++) {
          chi_prompt[p][a][0][gout] /= chi_sum;
        }

        // And extend to the remaining incoming groups
        for (int gin = 1; gin < energy_groups; gin++) {
          chi_prompt[p][a][gin] = chi_prompt[p][a][0];
        }

        // Finally set chi-delayed equal to chi-prompt
        // Set chi-delayed to chi-prompt
        for(int gin = 0; gin < energy_groups; gin++) {
          for (int gout = 0; gout < energy_groups; gout++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              chi_delayed[p][a][gin][gout][dg] =
                   chi_prompt[p][a][gin][gout];
            }
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
    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // nu-fission is a 3-d array
      read_nd_vector(xsdata_grp, "nu-fission", prompt_nu_fission);

      // set delayed-nu-fission and correct prompt-nu-fission with beta
      for (int p = 0; p < n_pol; p++) {
        for (int a = 0; a < n_azi; a++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              delayed_nu_fission[p][a][gin][dg] =
                   temp_beta[p][a][gin][dg] *
                   prompt_nu_fission[p][a][gin];
            }

            // Correct the prompt-nu-fission using the delayed neutron fraction
            if (delayed_groups > 0) {
              double beta_sum = std::accumulate(temp_beta[p][a][gin].begin(),
                                                temp_beta[p][a][gin].end(), 0.);
              prompt_nu_fission[p][a][gin] *= (1. - beta_sum);
            }
          }
        }
      }

    } else if (ndims == 4) {
      // nu-fission is a matrix
      read_nd_vector(xsdata_grp, "nu_fission", chi_prompt);

      // Normalize the chi info so the CDF is 1.
      for (int p = 0; p < n_pol; p++) {
        for (int a = 0; a < n_azi; a++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            double chi_sum = std::accumulate(chi_prompt[p][a][gin].begin(),
                                             chi_prompt[p][a][gin].end(), 0.);
            if (chi_sum >= 0.) {
              for (int gout = 0; gout < energy_groups; gout++) {
                chi_prompt[p][a][gin][gout] /= chi_sum;
              }
            } else {
              fatal_error("Encountered chi for a group that is <= 0!");
            }
          }

          // set chi-delayed to chi-prompt
          for (int gin = 0; gin < energy_groups; gin++) {
            for (int gout = 0; gout < energy_groups; gout++) {
              for (int dg = 0; dg < delayed_groups; dg++) {
                chi_delayed[p][a][gin][gout][dg] =
                     chi_prompt[p][a][gin][gout];
              }
            }
          }

          // Set the vector nu-fission from the matrix nu-fission
          for (int gin = 0; gin < energy_groups; gin++) {
            double sum = std::accumulate(chi_prompt[p][a][gin].begin(),
                                         chi_prompt[p][a][gin].end(), 0.);
            prompt_nu_fission[p][a][gin] = sum;
          }

          // Set the delayed-nu-fission and correct prompt-nu-fission with beta
          for (int gin = 0; gin < energy_groups; gin++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              delayed_nu_fission[p][a][gin][dg] =
                   temp_beta[p][a][gin][dg] *
                   prompt_nu_fission[p][a][gin];
            }

            // Correct prompt-nu-fission using the delayed neutron fraction
            if (delayed_groups > 0) {
              double beta_sum = std::accumulate(temp_beta[p][a][gin].begin(),
                                                temp_beta[p][a][gin].end(), 0.);
              prompt_nu_fission[p][a][gin] *= (1. - beta_sum);
            }
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
    double_3dvec temp_arr = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups)));
    read_nd_vector(xsdata_grp, "chi-prompt", temp_arr);

    for (int a = 0; a < n_azi; a++) {
      for (int p = 0; p < n_pol; p++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          for (int gout = 0; gout < energy_groups; gout++) {
            chi_prompt[p][a][gin][gout] = temp_arr[p][a][gout];
          }

          // Normalize chi so its CDF goes to 1
          double chi_sum = std::accumulate(chi_prompt[p][a][gin].begin(),
                                           chi_prompt[p][a][gin].end(), 0.);
          if (chi_sum >= 0.) {
            for (int gout = 0; gout < energy_groups; gout++) {
              chi_prompt[p][a][gin][gout] /= chi_sum;
            }
          } else {
            fatal_error("Encountered chi-prompt for a group that is <= 0.!");
          }
        }
      }
    }
  }

  // If chi-delayed is provided, set chi-delayed
  if (object_exists(xsdata_grp, "chi-delayed")) {
    hid_t xsdata = open_dataset(xsdata_grp, "chi-delayed");
    int ndims = dataset_ndims(xsdata);
    if (is_isotropic) ndims += 2;
    close_dataset(xsdata);

    if (ndims == 3) {
      // chi-delayed is a [in group] vector
      double_3dvec temp_arr = double_3dvec(n_pol, double_2dvec(n_azi,
           double_1dvec(energy_groups)));
      read_nd_vector(xsdata_grp, "chi-delayed", temp_arr);

      for (int a = 0; a < n_azi; a++) {
        for (int p = 0; p < n_pol; p++) {
          // normalize the chi CDF to 1
          double chi_sum = std::accumulate(temp_arr[p][a].begin(),
                                           temp_arr[p][a].end(), 0.);
          if (chi_sum <= 0.) {
            fatal_error("Encountered chi-delayed for a group that is <= 0!");
          }

          // set chi-delayed
          for (int gin = 0; gin < energy_groups; gin++) {
            for (int gout = 0; gout < energy_groups; gout++) {
              for (int dg = 0; dg < delayed_groups; dg++) {
                chi_delayed[p][a][gin][gout][dg] =
                     temp_arr[p][a][gout] / chi_sum;
              }
            }
          }
        }
      }
    } else if (ndims == 4) {
      // chi_delayed is a matrix
      read_nd_vector(xsdata_grp, "chi-delayed", chi_delayed);

      // Normalize the chi info so the CDF is 1.
      for (int a = 0; a < n_azi; a++) {
        for (int p = 0; p < n_pol; p++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            for (int gin = 0; gin < energy_groups; gin++) {
              double chi_sum = 0.;
              for (int gout = 0; gout < energy_groups; gout++) {
                chi_sum += chi_delayed[p][a][gin][gout][dg];
              }

              if (chi_sum > 0.) {
                for (int gout = 0; gout < energy_groups; gout++) {
                  chi_delayed[p][a][gin][gout][dg] /= chi_sum;
                }
              } else {
                fatal_error("Encountered chi-delayed for a group that is <= 0!");
              }
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
    if (is_isotropic) ndims += 2;
    close_dataset(xsdata);

    if (ndims == 3) {
      // prompt-nu-fission is a [in group] vector
      read_nd_vector(xsdata_grp, "prompt-nu-fission",
                     prompt_nu_fission);
    } else if (ndims == 4) {
      // prompt nu fission is a matrix,
      // so set prompt_nu_fiss & chi_prompt
      double_4dvec temp_arr = double_4dvec(n_pol, double_3dvec(n_azi,
           double_2dvec(energy_groups, double_1dvec(energy_groups))));
      read_nd_vector(xsdata_grp, "prompt-nu-fission", temp_arr);

      // The prompt_nu_fission vector from the matrix form
      for (int a = 0; a < n_azi; a++) {
        for (int p = 0; p < n_pol; p++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            double prompt_sum = std::accumulate(temp_arr[p][a][gin].begin(),
                                                temp_arr[p][a][gin].end(), 0.);
            prompt_nu_fission[p][a][gin] = prompt_sum;
          }

          // The chi_prompt data is just the normalized fission matrix
          for (int gin= 0; gin < energy_groups; gin++) {
            if (prompt_nu_fission[p][a][gin] > 0.) {
              for (int gout = 0; gout < energy_groups; gout++) {
                chi_prompt[p][a][gin][gout] =
                     temp_arr[p][a][gin][gout] /
                     prompt_nu_fission[p][a][gin];
              }
            } else {
              fatal_error("Encountered chi-prompt for a group that is <= 0!");
            }
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
    if (is_isotropic) ndims += 2;

    if (ndims == 3) {
      // delayed-nu-fission is a [in group] vector
      if (temp_beta[0][0][0][0] == 0.) {
        fatal_error("cannot set delayed-nu-fission with a 1D array if "
                    "beta is not provided");
      }
      double_3dvec temp_arr = double_3dvec(n_pol, double_2dvec(n_azi,
         double_1dvec(energy_groups)));
      read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_arr);

      for (int p = 0; p < n_pol; p++) {
        for (int a = 0; a < n_azi; a++) {
          for (int gin = 0; gin < energy_groups; gin++) {
            for (int dg = 0; dg < delayed_groups; dg++) {
              // Set delayed-nu-fission using beta
              delayed_nu_fission[p][a][gin][dg] =
                   temp_beta[p][a][gin][dg] * temp_arr[p][a][gin];
            }
          }
        }
      }

    } else if (ndims == 4) {
      read_nd_vector(xsdata_grp, "delayed-nu-fission",
                    delayed_nu_fission);

    } else if (ndims == 5) {
      // This will contain delayed-nu-fision and chi-delayed data
      double_5dvec temp_arr = double_5dvec(n_pol, double_4dvec(n_azi,
         double_3dvec(energy_groups, double_2dvec(energy_groups,
         double_1dvec(delayed_groups)))));
      read_nd_vector(xsdata_grp, "delayed-nu-fission", temp_arr);

      // Set the 4D delayed-nu-fission matrix and 5D chi-delayed matrix
      // from the 5D delayed-nu-fission matrix
      for (int p = 0; p < n_pol; p++) {
        for (int a = 0; a < n_azi; a++) {
          for (int dg = 0; dg < delayed_groups; dg++) {
            for (int gin = 0; gin < energy_groups; gin++) {
              double gout_sum = 0.;
              for (int gout = 0; gout < energy_groups; gout++) {
                gout_sum += temp_arr[p][a][gin][gout][dg];
                chi_delayed[p][a][gin][gout][dg] =
                     temp_arr[p][a][gin][gout][dg];
              }
              delayed_nu_fission[p][a][gin][dg] = gout_sum;
              // Normalize chi-delayed
              if (gout_sum > 0.) {
                for (int gout = 0; gout < energy_groups; gout++) {
                  chi_delayed[p][a][gin][gout][dg] /= gout_sum;
                }
              } else {
                fatal_error("Encountered chi-delayed for a group that is <= 0!");
              }
            }
          }
        }
      }

    } else {
      fatal_error("prompt-nu-fission must be provided as a 3D, 4D, or 5D "
                  "array!");
    }
    close_dataset(xsdata);
  }

  // Get the fission and kappa_fission data xs
  read_nd_vector(xsdata_grp, "fission", fission);
  read_nd_vector(xsdata_grp, "kappa-fission", kappa_fission);
}


void XsData::_scatter_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
     int energy_groups, int scatter_format, int final_scatter_format,
     int order_data, int max_order, int legendre_to_tabular_points)
{
  if (!object_exists(xsdata_grp, "scatter_data")) {
    fatal_error("Must provide scatter_data group!");
  }
  hid_t scatt_grp = open_group(xsdata_grp, "scatter_data");

  // Get the outgoing group boundary indices
  int_3dvec gmin = int_3dvec(n_pol, int_2dvec(n_azi,
       int_1dvec(energy_groups)));
  read_nd_vector(scatt_grp, "g_min", gmin, true);
  int_3dvec gmax = int_3dvec(n_pol, int_2dvec(n_azi,
       int_1dvec(energy_groups)));
  read_nd_vector(scatt_grp, "g_max", gmax, true);

  // Make gmin and gmax start from 0 vice 1 as they do in the library
  for (int p = 0; p < n_pol; p++) {
    for (int a = 0; a < n_azi; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        gmin[p][a][gin] -= 1;
        gmax[p][a][gin] -= 1;
      }
    }
  }

  // Now use this info to find the length of a vector to hold the flattened
  // data.
  int length = 0;
  for (int p = 0; p < n_pol; p++) {
    for (int a = 0; a < n_azi; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        length += order_data * (gmax[p][a][gin] - gmin[p][a][gin] + 1);
      }
    }
  }
  double_1dvec temp_arr = double_1dvec(length);
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
  double_5dvec input_scatt =
       double_5dvec(n_pol, double_4dvec(n_azi, double_3dvec(energy_groups)));

  int temp_idx = 0;
  for (int p = 0; p < n_pol; p++) {
    for (int a = 0; a < n_azi; a++) {
      for (int gin = 0; gin < energy_groups; gin++) {
        input_scatt[p][a][gin].resize(gmax[p][a][gin] - gmin[p][a][gin] + 1);
        for (int i_gout = 0; i_gout < input_scatt[p][a][gin].size(); i_gout++) {
          input_scatt[p][a][gin][i_gout].resize(order_dim);
          for (int l = 0; l < order_dim; l++) {
            input_scatt[p][a][gin][i_gout][l] = temp_arr[temp_idx++];
          }
          // Adjust index for the orders we didnt take
          temp_idx += (order_data - order_dim);
        }
      }
    }
  }
  temp_arr.clear();

  // Get multiplication matrix
  double_4dvec temp_mult = double_4dvec(n_pol, double_3dvec(n_azi,
       double_2dvec(energy_groups)));
  if (object_exists(scatt_grp, "multiplicity_matrix")) {
    temp_arr.resize(length);
    read_nd_vector(scatt_grp, "multiplicity_matrix", temp_arr);

    // convert the flat temp_arr to a jagged array for passing to scatt data
    int temp_idx = 0;
    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          temp_mult[p][a][gin].resize(gmax[p][a][gin] - gmin[p][a][gin] + 1);
          for (int i_gout = 0; i_gout < temp_mult[p][a][gin].size(); i_gout++) {
            temp_mult[p][a][gin][i_gout] = temp_arr[temp_idx++];
          }
        }
      }
    }
  } else {
    // Use a default: multiplicities are 1.0.
    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        for (int gin = 0; gin < energy_groups; gin++) {
          temp_mult[p][a][gin].resize(gmax[p][a][gin] - gmin[p][a][gin] + 1);
          for (int i_gout = 0; i_gout < temp_mult[p][a][gin].size(); i_gout++) {
            temp_mult[p][a][gin][i_gout] = 1.;
          }
        }
      }
    }
  }
  close_group(scatt_grp);

  // Finally, convert the Legendre data to tabular, if needed
  if (scatter_format == ANGLE_LEGENDRE &&
      final_scatter_format == ANGLE_TABULAR) {
    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        ScattDataLegendre legendre_scatt;
        legendre_scatt.init(gmin[p][a], gmax[p][a], temp_mult[p][a],
                            input_scatt[p][a]);

        // Now create a tabular version of legendre_scatt
        convert_legendre_to_tabular(legendre_scatt,
             *static_cast<ScattDataTabular*>(scatter[p][a]),
             legendre_to_tabular_points);

        scatter_format = final_scatter_format;
      }
    }
  } else {
    // We are sticking with the current representation
    // Initialize the ScattData object with this data
    for (int p = 0; p < n_pol; p++) {
      for (int a = 0; a < n_azi; a++) {
        scatter[p][a]->init(gmin[p][a], gmax[p][a], temp_mult[p][a],
                                  input_scatt[p][a]);
      }
    }
  }
}


void XsData::combine(std::vector<XsData*> those_xs, double_1dvec& scalars)
{
  // Combine the non-scattering data
  for (int i = 0; i < those_xs.size(); i++) {
    XsData* that = those_xs[i];
    if (!equiv(*that)) fatal_error("Cannot combine the XsData objects!");
    double scalar = scalars[i];
    for (int p = 0; p < total.size(); p++) {
      for (int a = 0; a < total[p].size(); a++) {
        for (int gin = 0; gin < total[p][a].size(); gin++) {
          total[p][a][gin] += scalar * that->total[p][a][gin];
          absorption[p][a][gin] += scalar * that->absorption[p][a][gin];
          inverse_velocity[p][a][gin] +=
                 scalar * that->inverse_velocity[p][a][gin];

          prompt_nu_fission[p][a][gin] +=
               scalar * that->prompt_nu_fission[p][a][gin];
          kappa_fission[p][a][gin] +=
               scalar * that->kappa_fission[p][a][gin];
          fission[p][a][gin] +=
               scalar * that->fission[p][a][gin];

          for (int dg = 0; dg < delayed_nu_fission[p][a][gin].size(); dg++) {
            delayed_nu_fission[p][a][gin][dg] +=
                 scalar * that->delayed_nu_fission[p][a][gin][dg];
          }

          for (int gout = 0; gout < chi_prompt[p][a][gin].size(); gout++) {
            chi_prompt[p][a][gin][gout] +=
                 scalar * that->chi_prompt[p][a][gin][gout];

            for (int dg = 0; dg < chi_delayed[p][a][gin][gout].size(); dg++) {
              chi_delayed[p][a][gin][gout][dg] +=
                   scalar * that->chi_delayed[p][a][gin][gout][dg];
            }
          }
        }

        for (int dg = 0; dg < decay_rate[p][a].size(); dg++) {
          decay_rate[p][a][dg] += scalar * that->decay_rate[p][a][dg];
        }

        // Normalize chi
        for (int gin = 0; gin < chi_prompt[p][a].size(); gin++) {
          double norm = std::accumulate(chi_prompt[p][a][gin].begin(),
                                        chi_prompt[p][a][gin].end(), 0.);
          if (norm > 0.) {
            for (int gout = 0; gout < chi_prompt[p][a][gin].size(); gout++) {
              chi_prompt[p][a][gin][gout] /= norm;
            }
          }

          for (int dg = 0; dg < chi_delayed[p][a][gin][0].size(); dg++) {
            norm = 0.;
            for (int gout = 0; gout < chi_delayed[p][a][gin].size(); gout++) {
              norm += chi_delayed[p][a][gin][gout][dg];
            }
            if (norm > 0.) {
              for (int gout = 0; gout < chi_delayed[p][a][gin].size(); gout++) {
                chi_delayed[p][a][gin][gout][dg] /= norm;
              }
            }
          }
        }
      }
    }
  }

  // Allow the ScattData object to combine itself
  for (int p = 0; p < total.size(); p++) {
    for (int a = 0; a < total[p].size(); a++) {
      // Build vector of the scattering objects to incorporate
      std::vector<ScattData*> those_scatts(those_xs.size());
      for (int i = 0; i < those_xs.size(); i++) {
        those_scatts[i] = those_xs[i]->scatter[p][a];
      }

      // Now combine these guys
      scatter[p][a]->combine(those_scatts, scalars);
    }
  }
}


bool XsData::equiv(const XsData& that)
{
  bool match = false;
  // check n_pol (total.size()), n_azi (total[0].size()), and
  // groups (total[0][0].size())
  // This assumes correct initializatino of the remaining cross sections
  if ((total.size() == that.total.size()) &&
      (total[0].size() == that.total[0].size()) &&
      (total[0][0].size() == that.total[0][0].size())) {
    match = true;
  }
  return match;
}

} //namespace openmc
