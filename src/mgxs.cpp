#include "openmc/mgxs.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fmt/core.h>
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/mgxs_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"


namespace openmc {

//==============================================================================
// Mgxs base-class methods
//==============================================================================

void
Mgxs::init(const std::string& in_name, double in_awr,
     const std::vector<double>& in_kTs, bool in_fissionable,
     AngleDistributionType in_scatter_format, bool in_is_isotropic,
     const std::vector<double>& in_polar, const std::vector<double>& in_azimuthal)
{
  // Set the metadata
  name = in_name;
  awr = in_awr;
  //TODO: Remove adapt when in_KTs is an xtensor
  kTs = xt::adapt(in_kTs);
  fissionable = in_fissionable;
  scatter_format = in_scatter_format;
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
Mgxs::metadata_from_hdf5(hid_t xs_id, const std::vector<double>& temperature,
     std::vector<int>& temps_to_read, int& order_dim)
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
  size_t num_temps = get_num_datasets(kT_group);
  char** dset_names = new char*[num_temps];
  for (int i = 0; i < num_temps; i++) {
    dset_names[i] = new char[151];
  }
  get_datasets(kT_group, dset_names);
  std::vector<size_t> shape = {num_temps};
  xt::xarray<double> available_temps(shape);
  for (int i = 0; i < num_temps; i++) {
    read_double(kT_group, dset_names[i], &available_temps[i], true);

    // convert eV to Kelvin
    available_temps[i] /= K_BOLTZMANN;

    // Done with dset_names, so delete it
    delete[] dset_names[i];
  }
  delete[] dset_names;
  std::sort(available_temps.begin(), available_temps.end());

  // If only one temperature is available, lets just use nearest temperature
  // interpolation
  if ((num_temps == 1) && (settings::temperature_method == TemperatureMethod::INTERPOLATION)) {
    warning("Cross sections for " + strtrim(name) + " are only available " +
            "at one temperature.  Reverting to the nearest temperature " +
            "method.");
    settings::temperature_method = TemperatureMethod::NEAREST;
  }

  switch(settings::temperature_method) {
    case TemperatureMethod::NEAREST:
      // Determine actual temperatures to read
      for (const auto& T : temperature) {
        auto i_closest = xt::argmin(xt::abs(available_temps - T))[0];
        double temp_actual = available_temps[i_closest];
        if (std::fabs(temp_actual - T) < settings::temperature_tolerance) {
          if (std::find(temps_to_read.begin(), temps_to_read.end(), std::round(temp_actual))
              == temps_to_read.end()) {
            temps_to_read.push_back(std::round(temp_actual));
          }
        } else {
          fatal_error(fmt::format(
            "MGXS library does not contain cross sections for {} at or near {} K.",
            in_name, std::round(T)));
        }
      }
      break;

    case TemperatureMethod::INTERPOLATION:
      for (int i = 0; i < temperature.size(); i++) {
        for (int j = 0; j < num_temps; j++) {
          if (j == (num_temps - 1)) {
            fatal_error("MGXS Library does not contain cross sections for " +
                  in_name + " at temperatures that bound " +
                  std::to_string(std::round(temperature[i])));
          }
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
            break;
          }
        }

      }
  }
  std::sort(temps_to_read.begin(), temps_to_read.end());

  // Get the library's temperatures
  int n_temperature = temps_to_read.size();
  std::vector<double> in_kTs(n_temperature);
  for (int i = 0; i < n_temperature; i++) {
    std::string temp_str(std::to_string(temps_to_read[i]) + "K");

    //read exact temperature value
    read_double(kT_group, temp_str.c_str(), &in_kTs[i], true);
  }
  close_group(kT_group);

  // Load the remaining metadata
  AngleDistributionType in_scatter_format;
  if (attribute_exists(xs_id, "scatter_format")) {
    std::string temp_str(MAX_WORD_LEN, ' ');
    read_attr_string(xs_id, "scatter_format", MAX_WORD_LEN, &temp_str[0]);
    to_lower(strtrim(temp_str));
    if (temp_str.compare(0, 8, "legendre") == 0) {
      in_scatter_format = AngleDistributionType::LEGENDRE;
    } else if (temp_str.compare(0, 9, "histogram") == 0) {
      in_scatter_format = AngleDistributionType::HISTOGRAM;
    } else if (temp_str.compare(0, 7, "tabular") == 0) {
      in_scatter_format = AngleDistributionType::TABULAR;
    } else {
      fatal_error("Invalid scatter_format option!");
    }
  } else {
    in_scatter_format = AngleDistributionType::LEGENDRE;
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
  if (in_scatter_format == AngleDistributionType::LEGENDRE) {
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
  std::vector<double> in_polar(in_n_pol);
  double dangle = PI / in_n_pol;
  for (int p = 0; p < in_n_pol; p++) {
    in_polar[p]  = (p + 0.5) * dangle;
  }
  std::vector<double> in_azimuthal(in_n_azi);
  dangle = 2. * PI / in_n_azi;
  for (int a = 0; a < in_n_azi; a++) {
    in_azimuthal[a] = (a + 0.5) * dangle - PI;
  }

  // Finally use this data to initialize the MGXS Object
  init(in_name, in_awr, in_kTs, in_fissionable, in_scatter_format,
       in_is_isotropic, in_polar,
       in_azimuthal);
}

//==============================================================================

Mgxs::Mgxs(hid_t xs_id, const std::vector<double>& temperature,
    int num_group, int num_delay) :
  num_groups(num_group),
  num_delayed_groups(num_delay)
{
  // Call generic data gathering routine (will populate the metadata)
  int order_data;
  std::vector<int> temps_to_read;
  metadata_from_hdf5(xs_id, temperature, temps_to_read, order_data);

  // Set number of energy and delayed groups
  AngleDistributionType final_scatter_format = scatter_format;
  if (settings::legendre_to_tabular) {
    if (scatter_format == AngleDistributionType::LEGENDRE) final_scatter_format = AngleDistributionType::TABULAR;
  }

  // Load the more specific XsData information
  for (int t = 0; t < temps_to_read.size(); t++) {
    xs[t] = XsData(fissionable, final_scatter_format, n_pol, n_azi,
                   num_groups, num_delayed_groups);
    // Get the temperature as a string and then open the HDF5 group
    std::string temp_str = std::to_string(temps_to_read[t]) + "K";
    hid_t xsdata_grp = open_group(xs_id, temp_str.c_str());

    xs[t].from_hdf5(xsdata_grp, fissionable, scatter_format,
                    final_scatter_format, order_data, is_isotropic, n_pol, n_azi);
    close_group(xsdata_grp);

  } // end temperature loop

  // Make sure the scattering format is updated to the final case
  scatter_format = final_scatter_format;
}

//==============================================================================

Mgxs::Mgxs(const std::string& in_name, const std::vector<double>& mat_kTs,
     const std::vector<Mgxs*>& micros, const std::vector<double>& atom_densities,
     int num_group, int num_delay) :
  num_groups(num_group),
  num_delayed_groups(num_delay)
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
  AngleDistributionType in_scatter_format = micros[0]->scatter_format;
  bool in_is_isotropic = micros[0]->is_isotropic;
  std::vector<double> in_polar = micros[0]->polar;
  std::vector<double> in_azimuthal = micros[0]->azimuthal;

  init(in_name, in_awr, mat_kTs, in_fissionable, in_scatter_format,
       in_is_isotropic, in_polar, in_azimuthal);

  // Create the xs data for each temperature
  for (int t = 0; t < mat_kTs.size(); t++) {
    xs[t] = XsData(in_fissionable, in_scatter_format, in_polar.size(),
      in_azimuthal.size(), num_groups, num_delayed_groups);

    // Find the right temperature index to use
    double temp_desired = mat_kTs[t];

    // Create the list of temperature indices and interpolation factors for
    // each microscopic data at the material temperature
    std::vector<int> micro_t(micros.size(), 0);
    std::vector<double> micro_t_interp(micros.size(), 0.);
    for (int m = 0; m < micros.size(); m++) {
      switch(settings::temperature_method) {
      case TemperatureMethod::NEAREST:
        {
          micro_t[m] = xt::argmin(xt::abs(micros[m]->kTs - temp_desired))[0];
          auto temp_actual = micros[m]->kTs[micro_t[m]];

          if (std::abs(temp_actual - temp_desired) >= K_BOLTZMANN * settings::temperature_tolerance) {
            fatal_error(fmt::format(
              "MGXS Library does not contain cross section for {} at or near {} K.",
              name, std::round(temp_desired / K_BOLTZMANN)));
          }
        }
        break;
      case TemperatureMethod::INTERPOLATION:
        // Get a list of bounding temperatures for each actual temperature
        // present in the model
        for (int k = 0; k < micros[m]->kTs.shape()[0] - 1; k++) {
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
    if (settings::temperature_method == TemperatureMethod::NEAREST) num_interp_points = 1;
    std::vector<double> interp(micros.size());
    std::vector<int> temp_indices(micros.size());
    for (int interp_point = 0; interp_point < num_interp_points; interp_point++) {
      for (int m = 0; m < micros.size(); m++) {
        interp[m] = (1. - micro_t_interp[m]) * atom_densities[m];
        temp_indices[m] = micro_t[m] + interp_point;
        micro_t_interp[m] = 1. - micro_t_interp[m];
      }

      combine(micros, interp, temp_indices, t);
    } // end loop to sum all micros across the temperatures
  } // end temperature (t) loop
}

//==============================================================================

void
Mgxs::combine(const std::vector<Mgxs*>& micros, const std::vector<double>& scalars,
              const std::vector<int>& micro_ts, int this_t)
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
Mgxs::get_xs(MgxsType xstype, int gin, const int* gout, const double* mu,
  const int* dg)
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
  case MgxsType::TOTAL:
    val = xs_t->total(a, gin);
    break;
  case MgxsType::NU_FISSION:
    val = fissionable ? xs_t->nu_fission(a, gin) : 0.;
    break;
  case MgxsType::ABSORPTION:
    val = xs_t->absorption(a, gin);;
    break;
  case MgxsType::FISSION:
    val = fissionable ? xs_t->fission(a, gin) : 0.;
    break;
  case MgxsType::KAPPA_FISSION:
    val = fissionable ? xs_t->kappa_fission(a, gin) : 0.;
    break;
  case MgxsType::SCATTER:
  case MgxsType::SCATTER_MULT:
  case MgxsType::SCATTER_FMU_MULT:
  case MgxsType::SCATTER_FMU:
    val = xs_t->scatter[a]->get_xs(xstype, gin, gout, mu);
    break;
  case MgxsType::PROMPT_NU_FISSION:
    val = fissionable ? xs_t->prompt_nu_fission(a, gin) : 0.;
    break;
  case MgxsType::DELAYED_NU_FISSION:
    if (fissionable) {
      if (dg != nullptr) {
        val = xs_t->delayed_nu_fission(a, *dg, gin);
      } else {
        val = 0.;
        for (int d = 0; d < xs_t->delayed_nu_fission.shape()[1]; d++) {
          val += xs_t->delayed_nu_fission(a, d, gin);
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MgxsType::CHI_PROMPT:
    if (fissionable) {
      if (gout != nullptr) {
        val = xs_t->chi_prompt(a, gin, *gout);
      } else {
        // provide an outgoing group-wise sum
        val = 0.;
        for (int g = 0; g < xs_t->chi_prompt.shape()[2]; g++) {
            val += xs_t->chi_prompt(a, gin, g);
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MgxsType::CHI_DELAYED:
    if (fissionable) {
      if (gout != nullptr) {
        if (dg != nullptr) {
          val = xs_t->chi_delayed(a, *dg, gin, *gout);
        } else {
          val = xs_t->chi_delayed(a, 0, gin, *gout);
        }
      } else {
        if (dg != nullptr) {
          val = 0.;
          for (int g = 0; g < xs_t->delayed_nu_fission.shape()[2]; g++) {
            val += xs_t->delayed_nu_fission(a, *dg, gin, g);
          }
        } else {
          val = 0.;
          for (int g = 0; g < xs_t->delayed_nu_fission.shape()[2]; g++) {
            for (int d = 0; d < xs_t->delayed_nu_fission.shape()[3]; d++) {
             val += xs_t->delayed_nu_fission(a, d, gin, g);
            }
          }
        }
      }
    } else {
      val = 0.;
    }
    break;
  case MgxsType::INVERSE_VELOCITY:
    val = xs_t->inverse_velocity(a, gin);
    break;
  case MgxsType::DECAY_RATE:
    if (dg != nullptr) {
      val = xs_t->decay_rate(a, *dg);
    } else {
      val = xs_t->decay_rate(a, 0);
    }
    break;
  default:
    val = 0.;
  }
  return val;
}

//==============================================================================

void
Mgxs::sample_fission_energy(int gin, int& dg, int& gout, uint64_t* seed)
{
  // This method assumes that the temperature and angle indices are set
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  XsData* xs_t = &xs[cache[tid].t];
  double nu_fission = xs_t->nu_fission(cache[tid].a, gin);

  // Find the probability of having a prompt neutron
  double prob_prompt = xs_t->prompt_nu_fission(cache[tid].a, gin);

  // sample random numbers
  double xi_pd = prn(seed) * nu_fission;
  double xi_gout = prn(seed);

  // Select whether the neutron is prompt or delayed
  if (xi_pd <= prob_prompt) {
    // the neutron is prompt

    // set the delayed group for the particle to be -1, indicating prompt
    dg = -1;

    // sample the outgoing energy group
    double prob_gout = 0.;
    for (gout = 0; gout < num_groups; ++gout) {
      prob_gout += xs_t->chi_prompt(cache[tid].a, gin, gout);
      if (xi_gout < prob_gout) break;
    }

  } else {
    // the neutron is delayed

    // get the delayed group
    for (dg = 0; dg < num_delayed_groups; ++dg) {
      prob_prompt += xs_t->delayed_nu_fission(cache[tid].a, dg, gin);
      if (xi_pd < prob_prompt) break;
    }

    // adjust dg in case of round-off error
    dg = std::min(dg, num_delayed_groups - 1);

    // sample the outgoing energy group
    double prob_gout = 0.;
    for (gout = 0; gout < num_groups; ++gout) {
      prob_gout += xs_t->chi_delayed(cache[tid].a, dg, gin, gout);
      if (xi_gout < prob_gout) break;
    }
  }
}

//==============================================================================

void
Mgxs::sample_scatter(int gin, int& gout, double& mu, double& wgt, uint64_t* seed)
{
  // This method assumes that the temperature and angle indices are set
  // Sample the data
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  xs[cache[tid].t].scatter[cache[tid].a]->sample(gin, gout, mu, wgt, seed);
}

//==============================================================================

void
Mgxs::calculate_xs(Particle& p)
{
  // Set our indices
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  set_temperature_index(p.sqrtkT_);
  set_angle_index(p.u_local());
  XsData* xs_t = &xs[cache[tid].t];
  p.macro_xs_.total = xs_t->total(cache[tid].a, p.g_);
  p.macro_xs_.absorption = xs_t->absorption(cache[tid].a, p.g_);
  p.macro_xs_.nu_fission =
    fissionable ? xs_t->nu_fission(cache[tid].a, p.g_) : 0.;
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
    cache[tid].t = xt::argmin(xt::abs(kTs - sqrtkT * sqrtkT))[0];
    cache[tid].sqrtkT = sqrtkT;
  }
}

//==============================================================================

void
Mgxs::set_angle_index(Direction u)
{
  // See if we need to find the new index
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  if (!is_isotropic &&
      ((u.x != cache[tid].u) || (u.y != cache[tid].v) ||
       (u.z != cache[tid].w))) {
    // convert direction to polar and azimuthal angles
    double my_pol = std::acos(u.z);
    double my_azi = std::atan2(u.y, u.x);

    // Find the location, assuming equal-bin angles
    double delta_angle = PI / n_pol;
    int p = std::floor(my_pol / delta_angle);
    delta_angle = 2. * PI / n_azi;
    int a = std::floor((my_azi + PI) / delta_angle);

    cache[tid].a = n_azi * p + a;

    // store this direction as the last one used
    cache[tid].u = u.x;
    cache[tid].v = u.y;
    cache[tid].w = u.z;
  }
}

} // namespace openmc
