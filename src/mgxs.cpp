#include "openmc/mgxs.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <valarray>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/string_functions.h"


namespace openmc {

// Storage for the MGXS data
std::vector<Mgxs> nuclides_MG;
std::vector<Mgxs> macro_xs;


//==============================================================================
// Mgxs base-class methods
//==============================================================================

void
Mgxs::init(const std::string& in_name, double in_awr,
     const double_1dvec& in_kTs, bool in_fissionable, int in_scatter_format,
     int in_num_groups, int in_num_delayed_groups, bool in_is_isotropic,
     const double_1dvec& in_polar, const double_1dvec& in_azimuthal)
{
  // Set the metadata
  name = in_name;
  awr = in_awr;
  kTs = in_kTs;
  fissionable = in_fissionable;
  scatter_format = in_scatter_format;
  num_groups = in_num_groups;
  num_delayed_groups = in_num_delayed_groups;
  xs.resize(in_kTs.size());
  is_isotropic = in_is_isotropic;
  n_pol = in_polar.size();
  n_azi = in_azimuthal.size();
  polar = in_polar;
  azimuthal = in_azimuthal;

  // Set the cross section index cache
#ifdef _OPENMP
  int n_threads = omp_get_max_threads();
#else
  int n_threads = 1;
#endif
  cache.resize(n_threads);
  // std::vector.resize() will value-initialize the members of cache[:]
}

//==============================================================================

void
Mgxs::metadata_from_hdf5(hid_t xs_id, int in_num_groups,
     int in_num_delayed_groups, const double_1dvec& temperature,
     double tolerance, int_1dvec& temps_to_read, int& order_dim, int& method)
{
  // get name
  char char_name[MAX_WORD_LEN];
  get_name(xs_id, char_name);
  std::string in_name {char_name};
  // remove the leading '/'
  in_name = in_name.substr(1);

  // Get the AWR
  double in_awr;
  if (attribute_exists(xs_id, "atomic_weight_ratio")) {
    read_attr_double(xs_id, "atomic_weight_ratio", &in_awr);
  } else {
    in_awr = MACROSCOPIC_AWR;
  }

  // Determine the available temperatures
  hid_t kT_group = open_group(xs_id, "kTs");
  int num_temps = get_num_datasets(kT_group);
  char* dset_names[num_temps];
  for (int i = 0; i < num_temps; i++) {
    dset_names[i] = new char[151];
  }
  get_datasets(kT_group, dset_names);
  double_1dvec available_temps(num_temps);
  for (int i = 0; i < num_temps; i++) {
    read_double(kT_group, dset_names[i], &available_temps[i], true);

    // convert eV to Kelvin
    available_temps[i] /= K_BOLTZMANN;

    // Done with dset_names, so delete it
    delete[] dset_names[i];
  }
  std::sort(available_temps.begin(), available_temps.end());

  // If only one temperature is available, lets just use nearest temperature
  // interpolation
  if ((num_temps == 1) && (method == TEMPERATURE_INTERPOLATION)) {
    warning("Cross sections for " + strtrim(name) + " are only available " +
            "at one temperature.  Reverting to the nearest temperature " +
            "method.");
    method = TEMPERATURE_NEAREST;
  }

  switch(method) {
    case TEMPERATURE_NEAREST:
      // Find the minimum difference
      for (int i = 0; i < temperature.size(); i++) {
        std::valarray<double> temp_diff(available_temps.data(),
                                        available_temps.size());
        temp_diff = std::abs(temp_diff - temperature[i]);
        int i_closest = std::min_element(std::begin(temp_diff), std::end(temp_diff)) -
             std::begin(temp_diff);
        double temp_actual = available_temps[i_closest];

        if (std::abs(temp_actual - temperature[i]) < tolerance) {
          if (std::find(temps_to_read.begin(), temps_to_read.end(),
                        std::round(temp_actual)) == temps_to_read.end()) {
            temps_to_read.push_back(std::round(temp_actual));
          } else {
            fatal_error("MGXS Library does not contain cross section for " +
                        in_name + " at or near " +
                        std::to_string(std::round(temperature[i])) + " K.");
          }
        }
      }
      break;

    case TEMPERATURE_INTERPOLATION:
      for (int i = 0; i < temperature.size(); i++) {
        for (int j = 0; j < num_temps - 1; j++) {
          if ((available_temps[j] <= temperature[i]) &&
              (temperature[i] < available_temps[j + 1])) {
            if (std::find(temps_to_read.begin(),
                          temps_to_read.end(),
                          std::round(available_temps[j])) == temps_to_read.end()) {
              temps_to_read.push_back(std::round((int)available_temps[j]));
            }

            if (std::find(temps_to_read.begin(), temps_to_read.end(),
                          std::round(available_temps[j + 1])) == temps_to_read.end()) {
              temps_to_read.push_back(std::round((int) available_temps[j + 1]));
            }
            continue;
          }
        }

      fatal_error("MGXS Library does not contain cross sections for " +
                  in_name + " at temperatures that bound " +
                  std::to_string(std::round(temperature[i])));
      }
  }
  std::sort(temps_to_read.begin(), temps_to_read.end());

  // Get the library's temperatures
  int n_temperature = temps_to_read.size();
  double_1dvec in_kTs(n_temperature);
  for (int i = 0; i < n_temperature; i++) {
    std::string temp_str(std::to_string(temps_to_read[i]) + "K");

    //read exact temperature value
    read_double(kT_group, temp_str.c_str(), &in_kTs[i], true);
  }
  close_group(kT_group);

  // Load the remaining metadata
  int in_scatter_format;
  if (attribute_exists(xs_id, "scatter_format")) {
    std::string temp_str(MAX_WORD_LEN, ' ');
    read_attr_string(xs_id, "scatter_format", MAX_WORD_LEN, &temp_str[0]);
    to_lower(strtrim(temp_str));
    if (temp_str.compare(0, 8, "legendre") == 0) {
      in_scatter_format = ANGLE_LEGENDRE;
    } else if (temp_str.compare(0, 9, "histogram") == 0) {
      in_scatter_format = ANGLE_HISTOGRAM;
    } else if (temp_str.compare(0, 7, "tabular") == 0) {
      in_scatter_format = ANGLE_TABULAR;
    } else {
      fatal_error("Invalid scatter_format option!");
    }
  } else {
    in_scatter_format = ANGLE_LEGENDRE;
  }

  if (attribute_exists(xs_id, "scatter_shape")) {
    std::string temp_str(MAX_WORD_LEN, ' ');
    read_attr_string(xs_id, "scatter_shape", MAX_WORD_LEN, &temp_str[0]);
    to_lower(strtrim(temp_str));
    if (temp_str.compare(0, 14, "[g][g\'][order]") != 0) {
      fatal_error("Invalid scatter_shape option!");
    }
  }

  bool in_fissionable = false;
  if (attribute_exists(xs_id, "fissionable")) {
    int int_fiss;
    read_attr_int(xs_id, "fissionable", &int_fiss);
    in_fissionable = int_fiss;
  } else {
    fatal_error("Fissionable element must be set!");
  }

  // Get the library's value for the order
  if (attribute_exists(xs_id, "order")) {
    read_attr_int(xs_id, "order", &order_dim);
  } else {
    fatal_error("Order must be provided!");
  }

  // Store the dimensionality of the data in order_dim.
  // For Legendre data, we usually refer to it as Pn where n is the order.
  // However Pn has n+1 sets of points (since you need to count the P0
  // moment). Adjust for that. Histogram and Tabular formats dont need this
  // adjustment.
  if (in_scatter_format == ANGLE_LEGENDRE) {
    ++order_dim;
  }

  // Get the angular information
  int in_n_pol;
  int in_n_azi;
  bool in_is_isotropic = true;
  if (attribute_exists(xs_id, "representation")) {
    std::string temp_str(MAX_WORD_LEN, ' ');
    read_attr_string(xs_id, "representation", MAX_WORD_LEN, &temp_str[0]);
    to_lower(strtrim(temp_str));
    if (temp_str.compare(0, 5, "angle") == 0) {
      in_is_isotropic = false;
    } else if (temp_str.compare(0, 9, "isotropic") != 0) {
      fatal_error("Invalid Data Representation!");
    }
  }

  if (!in_is_isotropic) {
    if (attribute_exists(xs_id, "num_polar")) {
      read_attr_int(xs_id, "num_polar", &in_n_pol);
    } else {
      fatal_error("num_polar must be provided!");
    }
    if (attribute_exists(xs_id, "num_azimuthal")) {
      read_attr_int(xs_id, "num_azimuthal", &in_n_azi);
    } else {
      fatal_error("num_azimuthal must be provided!");
    }
  } else {
    in_n_pol = 1;
    in_n_azi = 1;
  }

  // Set the angular bins to use equally-spaced bins
  double_1dvec in_polar(in_n_pol);
  double dangle = PI / in_n_pol;
  for (int p = 0; p < in_n_pol; p++) {
    in_polar[p]  = (p + 0.5) * dangle;
  }
  double_1dvec in_azimuthal(in_n_azi);
  dangle = 2. * PI / in_n_azi;
  for (int a = 0; a < in_n_azi; a++) {
    in_azimuthal[a] = (a + 0.5) * dangle - PI;
  }

  // Finally use this data to initialize the MGXS Object
  init(in_name, in_awr, in_kTs, in_fissionable, in_scatter_format,
       in_num_groups, in_num_delayed_groups, in_is_isotropic, in_polar,
       in_azimuthal);
}

//==============================================================================

Mgxs::Mgxs(hid_t xs_id, int energy_groups, int delayed_groups,
     const double_1dvec& temperature, double tolerance, int max_order,
     bool legendre_to_tabular, int legendre_to_tabular_points, int& method)
{
  // Call generic data gathering routine (will populate the metadata)
  int order_data;
  int_1dvec temps_to_read;
  metadata_from_hdf5(xs_id, energy_groups, delayed_groups, temperature,
       tolerance, temps_to_read, order_data, method);

  // Set number of energy and delayed groups
  int final_scatter_format = scatter_format;
  if (legendre_to_tabular) {
    if (scatter_format == ANGLE_LEGENDRE) final_scatter_format = ANGLE_TABULAR;
  }

  // Load the more specific XsData information
  for (int t = 0; t < temps_to_read.size(); t++) {
    xs[t] = XsData(energy_groups, delayed_groups, fissionable,
                   final_scatter_format, n_pol, n_azi);
    // Get the temperature as a string and then open the HDF5 group
    std::string temp_str = std::to_string(temps_to_read[t]) + "K";
    hid_t xsdata_grp = open_group(xs_id, temp_str.c_str());

    xs[t].from_hdf5(xsdata_grp, fissionable, scatter_format,
                    final_scatter_format, order_data, max_order,
                    legendre_to_tabular_points, is_isotropic, n_pol, n_azi);
    close_group(xsdata_grp);

  } // end temperature loop

  // Make sure the scattering format is updated to the final case
  scatter_format = final_scatter_format;
}

//==============================================================================

Mgxs::Mgxs(const std::string& in_name, const double_1dvec& mat_kTs,
     const std::vector<Mgxs*>& micros, const double_1dvec& atom_densities,
     double tolerance, int& method)
{
  // Get the minimum data needed to initialize:
  // Dont need awr, but lets just initialize it anyways
  double in_awr = -1.;
  // start with the assumption it is not fissionable
  bool in_fissionable = false;
  for (int m = 0; m < micros.size(); m++) {
    if (micros[m]->fissionable) in_fissionable = true;
  }
  // Force all of the following data to be the same; these will be verified
  // to be true later
  int in_scatter_format = micros[0]->scatter_format;
  int in_num_groups = micros[0]->num_groups;
  int in_num_delayed_groups = micros[0]->num_delayed_groups;
  bool in_is_isotropic = micros[0]->is_isotropic;
  double_1dvec in_polar = micros[0]->polar;
  double_1dvec in_azimuthal = micros[0]->azimuthal;

  init(in_name, in_awr, mat_kTs, in_fissionable, in_scatter_format,
       in_num_groups, in_num_delayed_groups, in_is_isotropic, in_polar,
       in_azimuthal);

  // Create the xs data for each temperature
  for (int t = 0; t < mat_kTs.size(); t++) {
    xs[t] = XsData(in_num_groups, in_num_delayed_groups, in_fissionable,
       in_scatter_format, in_polar.size(), in_azimuthal.size());

    // Find the right temperature index to use
    double temp_desired = mat_kTs[t];

    // Create the list of temperature indices and interpolation factors for
    // each microscopic data at the material temperature
    int_1dvec micro_t(micros.size(), 0);
    double_1dvec micro_t_interp(micros.size(), 0.);
    for (int m = 0; m < micros.size(); m++) {
      switch(method) {
      case TEMPERATURE_NEAREST:
        {
          // Find the nearest temperature
          std::valarray<double> temp_diff(micros[m]->kTs.data(),
                                          micros[m]->kTs.size());
          temp_diff = std::abs(temp_diff - temp_desired);
          micro_t[m] = std::min_element(std::begin(temp_diff),
                                        std::end(temp_diff)) -
               std::begin(temp_diff);
          double temp_actual = micros[m]->kTs[micro_t[m]];

          if (std::abs(temp_actual - temp_desired) >= K_BOLTZMANN * tolerance) {
            fatal_error("MGXS Library does not contain cross section for " +
                        name + " at or near " +
                        std::to_string(std::round(temp_desired / K_BOLTZMANN))
                        + " K.");
          }
        }
        break;
      case TEMPERATURE_INTERPOLATION:
        // Get a list of bounding temperatures for each actual temperature
        // present in the model
        for (int k = 0; k < micros[m]->kTs.size() - 1; k++) {
          if ((micros[m]->kTs[k] <= temp_desired) &&
              (temp_desired < micros[m]->kTs[k + 1])) {
            micro_t[m] = k;
            if (k == 0) {
              micro_t_interp[m] = (temp_desired - micros[m]->kTs[k]) /
                   (micros[m]->kTs[k + 1] - micros[m]->kTs[k]);
            } else {
              micro_t_interp[m] = 1.;
            }
          }
        }
      } // end switch
    } // end microscopic temperature loop

    // We are about to loop through each of the microscopic objects
    // and incorporate the contribution of each microscopic data at
    // one of the two temperature interpolants to this macroscopic quantity.
    // If we are doing nearest temperature interpolation, then we don't need
    // to do the 2nd temperature
    int num_interp_points = 2;
    if (method == TEMPERATURE_NEAREST) num_interp_points = 1;
    for (int interp_point = 0; interp_point < num_interp_points; interp_point++) {
      double_1dvec interp(micros.size());
      double_1dvec temp_indices(micros.size());
      for (int m = 0; m < micros.size(); m++) {
        interp[m] = (1. - micro_t_interp[m]) * atom_densities[m];
        temp_indices[m] = micro_t[m] + interp_point;
      }

      combine(micros, interp, micro_t, t);
    } // end loop to sum all micros across the temperatures
  } // end temperature (t) loop
}

//==============================================================================

void
Mgxs::combine(const std::vector<Mgxs*>& micros, const double_1dvec& scalars,
              const int_1dvec& micro_ts, int this_t)
{
  // Build the vector of pointers to the xs objects within micros
  std::vector<XsData*> those_xs(micros.size());
  for (int i = 0; i < micros.size(); i++) {
    if (!xs[this_t].equiv(micros[i]->xs[micro_ts[i]])) {
      fatal_error("Cannot combine the Mgxs objects!");
    }
    those_xs[i] = &(micros[i]->xs[micro_ts[i]]);
  }

  xs[this_t].combine(those_xs, scalars);
}

//==============================================================================

double
Mgxs::get_xs(int xstype, int gin, int* gout, double* mu, int* dg)
{
  // This method assumes that the temperature and angle indices are set
#ifdef _OPENMP
  int tid = omp_get_thread_num();
  XsData* xs_t = &xs[cache[tid].t];
  int a = cache[tid].a;
#else
  XsData* xs_t = &xs[cache[0].t];
  int a = cache[0].a;
#endif
  double val;
  switch(xstype) {
  case MG_GET_XS_TOTAL:
    val = xs_t->total[a][gin];
    break;
  case MG_GET_XS_NU_FISSION:
    val = fissionable ? xs_t->nu_fission[a][gin] : 0.;
    break;
  case MG_GET_XS_ABSORPTION:
    val = xs_t->absorption[a][gin];
    break;
  case MG_GET_XS_FISSION:
    val = fissionable ? xs_t->fission[a][gin] : 0.;
    break;
  case MG_GET_XS_KAPPA_FISSION:
    val = fissionable ? xs_t->kappa_fission[a][gin] : 0.;
    break;
  case MG_GET_XS_SCATTER:
  case MG_GET_XS_SCATTER_MULT:
  case MG_GET_XS_SCATTER_FMU_MULT:
  case MG_GET_XS_SCATTER_FMU:
    val = xs_t->scatter[a]->get_xs(xstype, gin, gout, mu);
    break;
  case MG_GET_XS_PROMPT_NU_FISSION:
    val = fissionable ? xs_t->prompt_nu_fission[a][gin] : 0.;
    break;
  case MG_GET_XS_DELAYED_NU_FISSION:
    if (fissionable) {
      if (dg != nullptr) {
        val = xs_t->delayed_nu_fission[a][gin][*dg];
      } else {
        val = 0.;
        for (auto& num : xs_t->delayed_nu_fission[a][gin]) {
          val += num;
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MG_GET_XS_CHI_PROMPT:
    if (fissionable) {
      if (gout != nullptr) {
        val = xs_t->chi_prompt[a][gin][*gout];
      } else {
        // provide an outgoing group-wise sum
        val = 0.;
        for (auto& num : xs_t->chi_prompt[a][gin]) {
          val += num;
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MG_GET_XS_CHI_DELAYED:
    if (fissionable) {
      if (gout != nullptr) {
        if (dg != nullptr) {
          val = xs_t->chi_delayed[a][gin][*gout][*dg];
        } else {
          val = xs_t->chi_delayed[a][gin][*gout][0];
        }
      } else {
        if (dg != nullptr) {
          val = 0.;
          for (int i = 0; i < xs_t->chi_delayed[a][gin].size(); i++) {
            val += xs_t->chi_delayed[a][gin][i][*dg];
          }
        } else {
          val = 0.;
          for (int i = 0; i < xs_t->chi_delayed[a][gin].size(); i++) {
            for (auto& num : xs_t->chi_delayed[a][gin][i]) {
              val += num;
            }
          }
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MG_GET_XS_INVERSE_VELOCITY:
    val = xs_t->inverse_velocity[a][gin];
    break;
  case MG_GET_XS_DECAY_RATE:
    if (dg != nullptr) {
      val = xs_t->decay_rate[a][*dg + 1];
    } else {
      val = xs_t->decay_rate[a][0];
    }
    break;
  default:
    val = 0.;
  }
  return val;
}

//==============================================================================

void
Mgxs::sample_fission_energy(int gin, int& dg, int& gout)
{
  // This method assumes that the temperature and angle indices are set
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  XsData* xs_t = &xs[cache[tid].t];
  double nu_fission = xs_t->nu_fission[cache[tid].a][gin];

  // Find the probability of having a prompt neutron
  double prob_prompt =
       xs_t->prompt_nu_fission[cache[tid].a][gin];

  // sample random numbers
  double xi_pd = prn() * nu_fission;
  double xi_gout = prn();

  // Select whether the neutron is prompt or delayed
  if (xi_pd <= prob_prompt) {
    // the neutron is prompt

    // set the delayed group for the particle to be -1, indicating prompt
    dg = -1;

    // sample the outgoing energy group
    gout = 0;
    double prob_gout =
         xs_t->chi_prompt[cache[tid].a][gin][gout];
    while (prob_gout < xi_gout) {
      gout++;
      prob_gout += xs_t->chi_prompt[cache[tid].a][gin][gout];
    }

  } else {
    // the neutron is delayed

    // get the delayed group
    dg = 0;
    while (xi_pd >= prob_prompt) {
      dg++;
      prob_prompt +=
           xs_t->delayed_nu_fission[cache[tid].a][gin][dg];
    }

    // adjust dg in case of round-off error
    dg = std::min(dg, num_delayed_groups - 1);

    // sample the outgoing energy group
    gout = 0;
    double prob_gout =
         xs_t->chi_delayed[cache[tid].a][gin][gout][dg];
    while (prob_gout < xi_gout) {
      gout++;
      prob_gout +=
           xs_t->chi_delayed[cache[tid].a][gin][gout][dg];
    }
  }
}

//==============================================================================

void
Mgxs::sample_scatter(int gin, int& gout, double& mu, double& wgt)
{
  // This method assumes that the temperature and angle indices are set
  // Sample the data
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  xs[cache[tid].t].scatter[cache[tid].a]->sample(gin, gout, mu, wgt);
}

//==============================================================================

void
Mgxs::calculate_xs(int gin, double sqrtkT, const double uvw[3],
     double& total_xs, double& abs_xs, double& nu_fiss_xs)
{
  // Set our indices
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  set_temperature_index(sqrtkT);
  set_angle_index(uvw);
  XsData* xs_t = &xs[cache[tid].t];
  total_xs = xs_t->total[cache[tid].a][gin];
  abs_xs = xs_t->absorption[cache[tid].a][gin];

  nu_fiss_xs = fissionable ? xs_t->nu_fission[cache[tid].a][gin] : 0.;
}

//==============================================================================

bool
Mgxs::equiv(const Mgxs& that)
{
  return ((num_delayed_groups == that.num_delayed_groups) &&
       (num_groups == that.num_groups) &&
       (n_pol == that.n_pol) &&
       (n_azi == that.n_azi) &&
       (std::equal(polar.begin(), polar.end(), that.polar.begin())) &&
       (std::equal(azimuthal.begin(), azimuthal.end(), that.azimuthal.begin())) &&
       (scatter_format == that.scatter_format));
}

//==============================================================================

void
Mgxs::set_temperature_index(double sqrtkT)
{
  // See if we need to find the new index
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  if (sqrtkT != cache[tid].sqrtkT) {
    double kT = sqrtkT * sqrtkT;

    // initialize vector for storage of the differences
    std::valarray<double> temp_diff(kTs.data(), kTs.size());

    // Find the minimum difference of kT and kTs
    temp_diff = std::abs(temp_diff - kT);
    cache[tid].t = std::min_element(std::begin(temp_diff), std::end(temp_diff)) -
         std::begin(temp_diff);

    // store this temperature as the last one used
    cache[tid].sqrtkT = sqrtkT;
  }
}

//==============================================================================

void
Mgxs::set_angle_index(const double uvw[3])
{
  // See if we need to find the new index
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  if (!is_isotropic &&
      ((uvw[0] != cache[tid].u) || (uvw[1] != cache[tid].v) ||
       (uvw[2] != cache[tid].w))) {
    // convert uvw to polar and azimuthal angles
    double my_pol = std::acos(uvw[2]);
    double my_azi = std::atan2(uvw[1], uvw[0]);

    // Find the location, assuming equal-bin angles
    double delta_angle = PI / n_pol;
    int p = std::floor(my_pol / delta_angle);
    delta_angle = 2. * PI / n_azi;
    int a = std::floor((my_azi + PI) / delta_angle);

    cache[tid].a = n_azi * p + a;

    // store this direction as the last one used
    cache[tid].u = uvw[0];
    cache[tid].v = uvw[1];
    cache[tid].w = uvw[2];
  }
}

} // namespace openmc
