#include "mgxs.h"

namespace openmc {

//==============================================================================
// Mgxs base-class methods
//==============================================================================

void Mgxs::init(const std::string& in_name, const double in_awr,
     const double_1dvec& in_kTs, const bool in_fissionable,
     const int in_scatter_format, const int in_num_groups,
     const int in_num_delayed_groups, const double_1dvec& in_polar,
     const double_1dvec& in_azimuthal)
{
  name = in_name;
  awr = in_awr;
  kTs = in_kTs;
  fissionable = in_fissionable;
  scatter_format = in_scatter_format;
  num_groups = in_num_groups;
  num_delayed_groups = in_num_delayed_groups;
  xs.resize(in_kTs.size());
  polar = in_polar;
  azimuthal = in_azimuthal;
  n_pol = polar.size();
  n_azi = azimuthal.size();
}


void Mgxs::_metadata_from_hdf5(const hid_t xs_id, const int in_num_groups,
     const int in_num_delayed_groups, double_1dvec& temperature, int& method,
     const double tolerance, int_1dvec& temps_to_read, int& order_dim,
     bool& is_isotropic)
{
  // get name
  char char_name[MAX_WORD_LEN];
  get_name(xs_id, char_name);
  std::string in_name(char_name, std::strlen(char_name));
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
  }
  std::sort(available_temps.begin(), available_temps.end());

  // If only one temperature is available, lets just use nearest temperature
  // interpolation
  if ((num_temps == 1) && (method == TEMPERATURE_INTERPOLATION)) {
    warning("Cross sections for " + strtrim(name) + " are only available " +
            "at one temperature.  Reverying to the nearest temperature " +
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
  //TODO: do i even need this flag? - it should be easy to self-determine
  bool in_fissionable = false;
  if (attribute_exists(xs_id, "fissionable")) {
    int int_fiss;
    read_attr_int(xs_id, "fissionable", &int_fiss);
    in_fissionable = (bool)int_fiss;
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
    order_dim = order_dim + 1;
  }

  // Get the angular information
  is_isotropic = true;
  if (attribute_exists(xs_id, "representation")) {
    std::string temp_str(MAX_WORD_LEN, ' ');
    read_attr_string(xs_id, "representation", MAX_WORD_LEN, &temp_str[0]);
    to_lower(strtrim(temp_str));
    if (temp_str.compare(0, 5, "angle") == 0) {
      is_isotropic =  false;
    } else if (temp_str.compare(0, 9, "isotropic") != 0) {
      fatal_error("Invalid Data Representation!");
    }
  }

  if (!is_isotropic) {
    if (attribute_exists(xs_id, "num_polar")) {
      read_attr_int(xs_id, "num_polar", &n_pol);
    } else {
      fatal_error("num_polar must be provided!");
    }
    if (attribute_exists(xs_id, "num_azimuthal")) {
      read_attr_int(xs_id, "num_azimuthal", &n_azi);
    } else {
      fatal_error("num_azimuthal must be provided!");
    }
  } else {
    n_pol = 1;
    n_azi = 1;
  }

  // Set the angular bins to use equally-spaced bins
  double_1dvec in_polar(n_pol);
  double dangle = PI / n_pol;
  for (int p = 0; p < n_pol; p++) {
    in_polar[p]  = (p + 0.5) * dangle;
  }
  double_1dvec in_azimuthal(n_azi);
  dangle = 2. * PI / n_azi;
  for (int a = 0; a < n_azi; a++) {
    in_azimuthal[a] = (a + 0.5) * dangle - PI;
  }

  // Finally use this data to initialize the MGXS Object
  init(in_name, in_awr, in_kTs, in_fissionable, in_scatter_format,
       in_num_groups, in_num_delayed_groups, in_polar, in_azimuthal);
}


void Mgxs::from_hdf5(hid_t xs_id, int energy_groups, int delayed_groups,
                     double_1dvec& temperature, int& method, double tolerance,
                     int max_order, bool legendre_to_tabular,
                     int legendre_to_tabular_points)
{
  // Call generic data gathering routine (will populate the metadata)
  int order_data;
  int_1dvec temps_to_read;
  bool is_isotropic;
  _metadata_from_hdf5(xs_id, energy_groups, delayed_groups, temperature,
       method, tolerance, temps_to_read, order_data, is_isotropic);

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
                    legendre_to_tabular_points, is_isotropic);
    close_group(xsdata_grp);

  } // end temperature loop
}


void Mgxs::build_macro(const std::string& in_name, double_1dvec& mat_kTs,
                       std::vector<Mgxs>& micros, double_1dvec& atom_densities,
                       int& method, double tolerance)
{
  // Get the minimum data needed to initialize:
  // Dont need awr, but lets just initialize it anyways
  double in_awr = -1.;
  // start with the assumption it is not fissionable
  bool in_fissionable = false;
  for (int m = 0; m < micros.size(); m++) {
    if (micros[m].fissionable) in_fissionable = true;
  }
  // Force all of the following data to be the same; these will be verified
  // to be true later
  int in_scatter_format = micros[0].scatter_format;
  int in_num_groups = micros[0].num_groups;
  int in_num_delayed_groups = micros[0].num_delayed_groups;
  double_1dvec in_polar = micros[0].polar;
  double_1dvec in_azimuthal = micros[0].azimuthal;

  init(in_name, in_awr, mat_kTs, in_fissionable, in_scatter_format, in_num_groups,
       in_num_delayed_groups, in_polar, in_azimuthal);

  // Create the xs data for each temperature
  for (int t = 0; t < mat_kTs.size(); t++) {
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
          std::valarray<double> temp_diff(micros[m].kTs.data(),
                                          micros[m].kTs.size());
          temp_diff = std::abs(temp_diff - temp_desired);
          micro_t[m] = std::min_element(std::begin(temp_diff),
                                        std::end(temp_diff)) -
               std::begin(temp_diff);
          double temp_actual = micros[m].kTs[micro_t[m]];

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
        for (int k = 0; k < micros[m].kTs.size() - 1; k++) {
          if ((micros[m].kTs[k] <= temp_desired) &&
              (temp_desired < micros[m].kTs[k + 1])) {
            micro_t[m] = k;
            if (k == 0) {
              micro_t_interp[m] = (temp_desired - micros[m].kTs[k]) /
                   (micros[m].kTs[k + 1] - micros[m].kTs[k]);
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


void Mgxs::combine(std::vector<Mgxs>& micros, double_1dvec& scalars,
                   int_1dvec& micro_ts, int this_t)
{
  // Build the vector of pointers to the xs objects within micros
  std::vector<XsData*> those_xs(micros.size());
  for (int i = 0; i < micros.size(); i++) {
    if (!xs[this_t].equiv(micros[i].xs[micro_ts[i]])) {
      fatal_error("Cannot combine the Mgxs objects!");
    }
    those_xs[i] = &(micros[i].xs[micro_ts[i]]);
  }

  xs[this_t].combine(those_xs, scalars);
}


double Mgxs::get_xs(const char* xstype, int gin, int* gout, double* mu, int* dg)
{
  // This method assumes that the temperature and angle indices are set
  double val;
  if (std::strcmp(xstype, "total")) {
    val = xs[index_temp].total[index_pol][index_azi][gin];
  } else if (std::strcmp(xstype, "absorption")) {
    val = xs[index_temp].absorption[index_pol][index_azi][gin];
  } else if (std::strcmp(xstype, "inverse-velocity")) {
    val = xs[index_temp].inverse_velocity[index_pol][index_azi][gin];
  } else if (std::strcmp(xstype, "decay rate")) {
    if (dg != nullptr) {
      val = xs[index_temp].decay_rate[index_pol][index_azi][*dg + 1];
    } else {
      val = xs[index_temp].decay_rate[index_pol][index_azi][0];
    }
  } else if ((std::strcmp(xstype, "scatter")) ||
             (std::strcmp(xstype, "scatter/mult")) ||
             (std::strcmp(xstype, "scatter*f_mu/mult")) ||
             (std::strcmp(xstype, "scatter*f_mu"))) {
    val = xs[index_temp].scatter[index_pol]
                    [index_azi]->get_xs(xstype, gin, gout, mu);
  } else if (fissionable && std::strcmp(xstype, "fission")) {
    val = xs[index_temp].fission[index_pol][index_azi][gin];
  } else if (fissionable && std::strcmp(xstype, "kappa-fission")) {
    val = xs[index_temp].kappa_fission[index_pol][index_azi][gin];
  } else if (fissionable && std::strcmp(xstype, "prompt-nu-fission")) {
    val = xs[index_temp].prompt_nu_fission[index_pol][index_azi][gin];
  } else if (fissionable && std::strcmp(xstype, "delayed-nu-fission")) {
    if (dg != nullptr) {
      val = xs[index_temp].delayed_nu_fission[index_pol][index_azi][gin][*dg];
    } else {
      val = 0.;
      for (auto& num : xs[index_temp].delayed_nu_fission[index_pol]
                                            [index_azi][gin]) {
        val += num;
      }
    }
  } else if (fissionable && std::strcmp(xstype, "nu-fission")) {
    val = xs[index_temp].prompt_nu_fission[index_pol][index_azi][gin];
    for (auto& num : xs[index_temp].delayed_nu_fission[index_pol]
                                          [index_azi][gin]) {
      val += num;
    }
  } else if (fissionable && std::strcmp(xstype, "chi-prompt")) {
    if (gout != nullptr) {
      val = xs[index_temp].chi_prompt[index_pol][index_azi][gin][*gout];
    } else {
      // provide an outgoing group-wise sum
      val = 0.;
      for (auto& num : xs[index_temp].chi_prompt[index_pol][index_azi][gin]) {
        val += num;
      }
    }
  } else if (fissionable && std::strcmp(xstype, "chi-delayed")) {
    if (gout != nullptr) {
      if (dg != nullptr) {
        val = xs[index_temp].chi_delayed[index_pol][index_azi][gin][*gout][*dg];
      } else {
        val = xs[index_temp].chi_delayed[index_pol][index_azi][gin][*gout][0];
      }
    } else {
      if (dg != nullptr) {
        val = 0.;
        for (int i = 0; i < xs[index_temp].chi_delayed[index_pol]
                                          [index_azi][gin].size(); i++) {
          val += xs[index_temp].chi_delayed[index_pol][index_azi][gin][i][*dg];
        }
      } else {
        val = 0.;
        for (int i = 0; i < xs[index_temp].chi_delayed[index_pol]
                                          [index_azi][gin].size(); i++) {
          for (auto& num : xs[index_temp].chi_delayed[index_pol]
                                         [index_azi][gin][i]) {
            val += num;
          }
        }
      }
    }
  } else {
    val = 0.;
  }
  return val;
}


void Mgxs::sample_fission_energy(int gin, double nu_fission, int& dg, int& gout)
{
  // This method assumes that the temperature and angle indices are set
  // Find the probability of having a prompt neutron
  double prob_prompt =
       xs[index_temp].prompt_nu_fission[index_pol][index_azi][gin] /
       nu_fission;

  // sample random numbers
  double xi_pd = prn();
  double xi_gout = prn();

  // Select whether the neutron is prompt or delayed
  if (xi_pd <= prob_prompt) {
    // the neutron is prompt

    // set the delayed group for the particle to be 0, indicating prompt
    dg = 0;

    // sample the outgoing energy group
    gout = 0;
    double prob_gout =
         xs[index_temp].chi_prompt[index_pol][index_azi][gin][gout];
    while (prob_gout < xi_gout) {
      gout++;
      prob_gout += xs[index_temp].chi_prompt[index_pol][index_azi][gin][gout];
    }

  } else {
    // the neutron is delayed

    // get the delayed group
    dg = 0;
    while (xi_pd >= prob_prompt) {
      dg++;
      prob_prompt +=
           xs[index_temp].delayed_nu_fission[index_pol][index_azi][gin][dg] /
           nu_fission;
    }

    // adjust dg in case of round-off error
    dg = std::min(dg, num_delayed_groups);

    // sample the outgoing energy group
    gout = 0;
    double prob_gout =
         xs[index_temp].chi_delayed[index_pol][index_azi][gin][gout][dg];
    while (prob_gout < xi_gout) {
      gout++;
      prob_gout +=
           xs[index_temp].chi_delayed[index_pol][index_azi][gin][gout][dg];
    }
  }
}


void Mgxs::sample_scatter(dir_arr& uvw, int gin, int& gout, double& mu,
                          double& wgt)
{
  // This method assumes that the temperature and angle indices are set
  // Sample the data
  xs[index_temp].scatter[index_pol][index_azi]->sample(gin, gout, mu, wgt);
}


void Mgxs::calculate_xs(int gin, double sqrtkT, dir_arr& uvw, double& total_xs,
                        double& abs_xs, double& nu_fiss_xs)
{
  // Set our indices
  set_temperature_index(sqrtkT);
  set_angle_index(uvw);
  total_xs = xs[index_temp].total[index_pol][index_azi][gin];
  abs_xs = xs[index_temp].absorption[index_pol][index_azi][gin];

  // nu-fission is made up of the prompt and all the delayed nu_fission data
  nu_fiss_xs = xs[index_temp].prompt_nu_fission[index_pol][index_azi][gin];
  for (auto& val : xs[index_temp].delayed_nu_fission[index_pol][index_azi][gin]) {
    nu_fiss_xs += val;
  }
}


bool Mgxs::equiv(const Mgxs& that)
{
  bool match = false;

  if ((num_delayed_groups == that.num_delayed_groups) &&
      (num_groups == that.num_groups) &&
      (n_pol == that.n_pol) &&
      (n_azi == that.n_azi) &&
      (std::equal(polar.begin(), polar.end(), that.polar.begin())) &&
      (std::equal(azimuthal.begin(), azimuthal.end(), that.azimuthal.begin())) &&
      (scatter_format == that.scatter_format)) {
    match = true;
  }
  return match;
}


inline void Mgxs::set_temperature_index(double sqrtkT)
{
  // See if we need to find the new index
  if (sqrtkT != last_sqrtkT) {
    double kT = sqrtkT * sqrtkT;

    // initialize vector for storage of the differences
    std::valarray<double> temp_diff(kTs.data(), kTs.size());

    // Find the minimum difference of kT and kTs
    temp_diff = std::abs(temp_diff - kT);
    index_temp = std::min_element(std::begin(temp_diff), std::end(temp_diff)) -
         std::begin(temp_diff);

    // store this temperature as the last one used
    last_sqrtkT = sqrtkT;
  }
}


inline void Mgxs::set_angle_index(dir_arr& uvw)
{
  // See if we need to find the new index
  if (uvw != last_uvw) {
    // convert uvw to polar and azimuthal angles
    double my_pol = std::acos(uvw[2]);
    double my_azi = std::atan2(uvw[1], uvw[0]);

    // Find the location, assuming equal-bin angles
    double delta_angle = PI / n_pol;
    index_pol = std::floor(my_pol / delta_angle + 1.);
    delta_angle = PI / n_azi;
    index_azi = std::floor((my_azi + PI) / delta_angle + 1.);

    // store this direction as the last one used
    last_uvw = uvw;
  }
}

//==============================================================================
// Mgxs data loading methods
//==============================================================================

void add_mgxs(hid_t file_id, char* name, int energy_groups,
     int delayed_groups, int n_temps, double temps[], int& method,
     double tolerance, int max_order, bool legendre_to_tabular,
     int legendre_to_tabular_points)
{
  //!! mgxs_data.F90 will be modified to just create the list of names
  //!! in the order needed
  // Convert temps to a vector for the from_hdf5 function
  double_1dvec temperature;
  temperature.assign(temps, temps + n_temps);

  // TODO: C++ replacement for write_message
  // write_message("Loading " + std::string(names[i]) + " data...", 6);

  // Check to make sure cross section set exists in the library
  hid_t xs_grp;
  if (object_exists(file_id, name)) {
    xs_grp = open_group(file_id, name);
  } else {
    fatal_error("Data for " + std::string(name) + " does not exist in "
                + "provided MGXS Library");
  }

  Mgxs mg;
  mg.from_hdf5(xs_grp, energy_groups, delayed_groups,
       temperature, method, tolerance, max_order, legendre_to_tabular,
       legendre_to_tabular_points);

  nuclides_MG.push_back(mg);
}


bool query_fissionable(const int n_nuclides, const int i_nuclides[])
{
  bool result = false;
  for (int i = 0; i < n_nuclides; i++) {
    if (nuclides_MG[i_nuclides[i]].fissionable) result = true;
  }
  return result;
}


void create_macro_xs(int n_materials, double_2dvec& mat_kTs,
                     std::vector<std::string>& mat_names,
                     double_1dvec& atom_densities, int& method,
                     double tolerance)
{
  // TODO mat_kTs needs to be converted from Fortran
  // it is currently an array of type(VectorReal), a wrapper should convert to
  // the vector.
  macro_xs.resize(n_materials);

  for (int m = 0; m < n_materials; m++)
  {
    if (mat_kTs[m].size() > 0) {
      macro_xs[m].build_macro(mat_names[m], mat_kTs[m], nuclides_MG,
                              atom_densities, method, tolerance);
    }
  }
}


} // namespace openmc