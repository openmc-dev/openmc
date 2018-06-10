#include "scattdata.h"

namespace openmc {

//==============================================================================
// ScattData base-class methods
//==============================================================================

void ScattData::generic_init(int order, int_1dvec in_gmin,
     int_1dvec in_gmax, double_2dvec in_energy, double_2dvec in_mult)
{
  int groups = in_energy.size();

  gmin = in_gmin;
  gmax = in_gmax;
  energy.resize(groups);
  mult.resize(groups);
  dist.resize(groups);

  for (int gin = 0; gin < groups; gin++) {
    // Make sure the energy is normalized
    double norm = std::accumulate(in_energy[gin].begin(),
                                  in_energy[gin].end(), 0.);

    if (norm != 0.) {
      for (auto& n : in_energy[gin]) n /= norm;
    }

    // Store the inputted data
    energy[gin] = in_energy[gin];
    mult[gin] = in_mult[gin];

    // Initialize the distribution data
    dist[gin].resize(in_gmax[gin] - in_gmin[gin] + 1);
    for (auto& v : dist[gin]) {
      v.resize(order);
      for (auto& n : v) n = 0.;
    }
  }
}


void ScattData::sample_energy(int gin, int& gout, int& i_gout)
{
  // Sample the outgoing group
  double xi = prn();

  i_gout = 0;
  gout = gmin[gin];
  double prob = energy[gin][i_gout];
  while((prob < xi) && (gout < gmax[gin])) {
    gout++;
    i_gout++;
    prob += energy[gin][i_gout];
  }
}


double ScattData::get_xs(const char* xstype, int gin, int* gout, double* mu)
{
  // Set the outgoing group offset index as needed
  int i_gout = 0;
  if (gout != nullptr) {
    // short circuit the function if gout is from a zero portion of the
    // scattering matrix
    if ((*gout < gmin[gin]) || (*gout >= gmax[gin])) { // > gmax?
      return 0.;
    }
    i_gout = *gout - gmin[gin];
  }

  double val = 0.;
  if (std::strcmp(xstype, "scatter")) {
    if (gout != nullptr) {
      val = scattxs[gin] * energy[gin][i_gout];
    } else {
      val = scattxs[gin];
    }
  } else if (std::strcmp(xstype, "scatter/mult")) {
    if (gout != nullptr) {
      val = scattxs[gin] * energy[gin][i_gout] / mult[gin][i_gout];
    } else {
      val = scattxs[gin] / std::inner_product(mult[gin].begin(),
                                              mult[gin].end(),
                                              energy[gin].begin(), 0.0);
    }
  } else if (std::strcmp(xstype, "scatter*f_mu/mult")) {
    if ((gout != nullptr) && (mu != nullptr)) {
      val = scattxs[gin] * energy[gin][i_gout] * calc_f(gin, *gout, *mu);
    } else {
      // This is not an expected path (asking for f_mu without asking for a
      // group or mu is not useful
      fatal_error("Invalid call to get_xs");
    }
  } else if (std::strcmp(xstype, "scatter*f_mu")) {
    if ((gout != nullptr) && (mu != nullptr)) {
      val = scattxs[gin] * energy[gin][i_gout] * calc_f(gin, *gout, *mu) /
           mult[gin][i_gout];
    } else {
      // This is not an expected path (asking for f_mu without asking for a
      // group or mu is not useful
      fatal_error("Invalid call to get_xs");
    }
  }
  return val;
}


//==============================================================================
// ScattDataLegendre methods
//==============================================================================

void ScattDataLegendre::init(int_1dvec& in_gmin, int_1dvec& in_gmax,
                             double_2dvec& in_mult, double_3dvec& coeffs)
{
  int groups = coeffs.size();
  int order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Get the scattering cross section value by summing the un-normalized P0
  // coefficient in the variable matrix over all outgoing groups.
  scattxs.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = in_gmax[gin] - in_gmin[gin] + 1;
    scattxs[gin] = 0.;
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      scattxs[gin] = std::accumulate(matrix[gin][i_gout].begin(),
                                     matrix[gin][i_gout].end(),
                                     scattxs[gin]);
    }
  }

  // Build the energy transfer matrix from data in the variable matrix while
  // also normalizing the variable matrix itself
  // (forcing the CDF of f(mu=1) == 1)
  double_2dvec in_energy;
  in_energy.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = in_gmax[gin] - in_gmin[gin] + 1;
    in_energy[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      double norm = matrix[gin][i_gout][0];
      in_energy[gin][i_gout] = norm;
      if (norm != 0.) {
        for (auto& n : matrix[gin][i_gout]) n /= norm;
      }
    }
  }

  // Initialize the base class attributes
  ScattData::generic_init(order, in_gmin, in_gmax, in_energy, in_mult);

  // Set the distribution (sdata.dist) values and initialize max_val
  max_val.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = gmax[gin] - gmin[gin] + 1;
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      dist[gin][i_gout] = matrix[gin][i_gout];
    }
    max_val[gin].resize(num_groups);
    for (auto& n : max_val[gin]) n = 0.;
  }

  // Now update the maximum value
  update_max_val();
}


void ScattDataLegendre::update_max_val()
{
  int groups = max_val.size();
  // Step through the polynomial with fixed number of points to identify the
  // maximal value
  int Nmu = 1001;
  double dmu = 2. / (Nmu - 1);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = gmax[gin] - gmin[gin] + 1;
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      for (int imu = 0; imu < Nmu; imu++) {
        double mu;
        if (imu == 0) {
          mu = -1.;
        } else if (imu == (Nmu - 1)) {
          mu = 1.;
        } else {
          mu = -1. + (imu - 1) * dmu;
        }

        // Calculate probability
        double f = evaluate_legendre_c(dist[gin][i_gout].size() - 1,
                                       dist[gin][i_gout].data(), mu);

        // if this is a new maximum, store it
        if (f > max_val[gin][i_gout]) max_val[gin][i_gout] = f;
      } // end imu loop

      // Since we may not have caught the true max, add 10% margin
      max_val[gin][i_gout] *= 1.1;
    }
  }
}


double ScattDataLegendre::calc_f(int gin, int gout, double mu)
{
  // TODO: gout >= or gout >?
  double f;
  if ((gout < gmin[gin]) || (gout >= gmax[gin])) {
    f = 0.;
  } else {
      // TODO: size() -1 or just size?
    int i_gout = gout - gmin[gin]; //TODO: + 1?
    f = evaluate_legendre_c(dist[gin][i_gout].size() - 1,
                            dist[gin][i_gout].data(), mu);
  }
  return f;
}


void ScattDataLegendre::sample(int gin, int& gout, double& mu, double& wgt)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout);

  // Now we can sample mu using the scattering kernel using rejection
  // sampling from a rectangular bounding box
  double M = max_val[gin][i_gout];
  int samples = 0;

  while(true) {
    double mu = 2. * prn() - 1.;
    double f = calc_f(gin, gout, mu);
    if (f > 0.) {
      double u = prn() * M;
      if (u <= f) break;
    }
    samples++;
    if (samples > MAX_SAMPLE) {
        fatal_error("Maximum number of Legendre expansion samples reached");
    }
  };

  // Update the weight to reflect neutron multiplicity
  wgt *= mult[gin][i_gout];
}


void ScattDataLegendre::combine(std::vector<ScattData*>& those_scatts,
                                double_1dvec& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  int max_order = 0;
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataLegendre* that = dynamic_cast<ScattDataLegendre*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    int that_order = that->get_order();
    if (that_order > max_order) max_order = that_order;
  }
  max_order++;  // Add one since this is a Legendre

  // Get the groups as a shorthand
  int groups = dynamic_cast<ScattDataLegendre*>(those_scatts[0])->energy.size();

  // Now allocate and zero our storage spaces
  double_3dvec this_matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(max_order, 0.)));
  double_2dvec mult_numer(groups, double_1dvec(groups, 0.));
  double_2dvec mult_denom(groups, double_1dvec(groups, 0.));

  // Build the dense scattering and multiplicity matrices
  // Get the multiplicity_matrix
  // To combine from nuclidic data we need to use the final relationship
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // Developed as follows:
  // mult_{gg'} = nuScatt{g,g'} / Scatt{g,g'}
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) / sum(N_i*scatt_{i,g,g'})
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // nuscatt_{i,g,g'} can be reconstructed from the energy and scattxs member
  // variables
  for (int i = 0; i < those_scatts.size(); i++) {
    ScattDataLegendre* that = dynamic_cast<ScattDataLegendre*>(those_scatts[i]);

    // Build the dense matrix for that object
    double_3dvec that_matrix = that->get_matrix(max_order);

    // Now add that to this for the scattering and multiplicity
    for (int gin = 0; gin < groups; gin++) {
      // Only spend time adding that's gmin to gmax data since the rest will
      // be zeros
      int i_gout = 0;
      for (int gout = that->gmin[gin]; gout <= that->gmax[gin]; gout++) {
        // Do the scattering matrix
        for (int l = 0; l < max_order; l++) {
          this_matrix[gin][gout][l] += scalars[i] * that_matrix[gin][gout][l];
        }

        // Incorporate that's contribution to the multiplicity matrix data
        double nuscatt = that->scattxs[gin] * that->energy[gin][i_gout];
        mult_numer[gin][gout] += scalars[i] * nuscatt;
        if (that->mult[gin][i_gout] > 0.) {
          mult_denom[gin][gout] += scalars[i] * nuscatt / that->mult[gin][i_gout];
        } else {
          mult_denom[gin][gout] += scalars[i];
        }
        i_gout++;
      }
    }
  }

  // Combine mult_numer and mult_denom into the combined multiplicity matrix
  double_2dvec this_mult(groups, double_1dvec(groups, 1.));
  for (int gin = 0; gin < groups; gin++) {
    for (int gout = 0; gout < groups; gout++) {
      if (mult_denom[gin][gout] > 0.) {
        this_mult[gin][gout] = mult_numer[gin][gout] / mult_denom[gin][gout];
      }
    }
  }
  mult_numer.clear();
  mult_denom.clear();

  // We have the data, now we need to convert to a jagged array and then use
  // the initialize function to store it on the object.
  int_1dvec in_gmin(groups);
  int_1dvec in_gmax(groups);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);
  for (int gin = 0; gin < groups; gin++) {
    // Find the minimum and maximum group boundaries
    int gmin_;
    for (gmin_ = 0; gmin_ < groups; gmin_++) {
      bool non_zero = std::all_of(this_matrix[gin][gmin_].begin(),
                                  this_matrix[gin][gmin_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }
    int gmax_;
    for (gmax_ = groups - 1; gmax_ >= 0; gmax_--) {
      bool non_zero = std::all_of(this_matrix[gin][gmax_].begin(),
                                  this_matrix[gin][gmax_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }

    // treat the case of all values being 0
    if (gmin_ > gmax_) {
      gmin_ = gin;
      gmax_ = gin;
    }

    // Store the group bounds
    in_gmin[gin] = gmin_;
    in_gmax[gin] = gmax_;

    // Store the data in the compressed format
    sparse_scatter[gin].resize(gmax_ - gmin_ + 1);
    sparse_mult[gin].resize(gmax_ - gmin_ + 1);
    int i_gout = 0;
    for (int gout = gmin_; gout <= gmax_; gout++) {
      sparse_scatter[gin][i_gout] = this_matrix[gin][gout];
      sparse_mult[gin][i_gout] = this_mult[gin][gout];
    }
  }

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}


double_3dvec ScattDataLegendre::get_matrix(int max_order)
{
  // Get the sizes and initialize the data to 0
  int groups = energy.size();
  int order_dim = max_order + 1;
  double_3dvec matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(order_dim, 0.)));

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix[gin][gout][l] = scattxs[gin] * energy[gin][i_gout] *
             dist[gin][i_gout][l];
      }
    }
  }
  return matrix;
}

//==============================================================================
// ScattDataHistogram methods
//==============================================================================

void ScattDataHistogram::init(int_1dvec& in_gmin, int_1dvec& in_gmax,
                            double_2dvec& in_mult, double_3dvec& coeffs)
{
  int groups = coeffs.size();
  int order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Get the scattering cross section value by summing the distribution
  // over all the histogram bins in angle and outgoing energy groups
  scattxs.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    scattxs[gin] = 0.;
    for (int i_gout = 0; i_gout < matrix[gin].size(); i_gout++) {
      scattxs[gin] += std::accumulate(matrix[gin][i_gout].begin(),
                                      matrix[gin][i_gout].end(), 0.);
    }
  }

  // Build the energy transfer matrix from data in the variable matrix
  double_2dvec in_energy;
  in_energy.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = in_gmax[gin] - in_gmin[gin] + 1;
    in_energy[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      double norm = std::accumulate(matrix[gin][i_gout].begin(),
                                    matrix[gin][i_gout].end(), 0.);
      in_energy[gin][i_gout] = norm;
      if (norm != 0.) {
        for (auto& n : matrix[gin][i_gout]) n /= norm;
      }
    }
  }

  // Initialize the base class attributes
  ScattData::generic_init(order, in_gmin, in_gmax, in_energy,
                                    in_mult);

  // Build the angular distributio mu values
  mu = double_1dvec(order);
  dmu = 2. / order;
  mu[0] = -1.;
  for (int imu = 1; imu < order; imu++) {
    mu[imu] = -1. + imu * dmu;
  }

  // Calculate f(mu) and integrate it so we can avoid rejection sampling
  fmu.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = gmax[gin] - gmin[gin] + 1;
    fmu[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      fmu[gin][i_gout].resize(order);
      // The variable matrix contains f(mu); so directly assign it
      fmu[gin][i_gout] = matrix[gin][i_gout];

      // Integrate the histogram
      dist[gin][i_gout][0] = dmu * matrix[gin][i_gout][0];
      for (int imu = 1; imu < order; imu++) {
        dist[gin][i_gout][imu] = dmu * matrix[gin][i_gout][imu] +
             dist[gin][i_gout][imu - 1];
      }

      // Now re-normalize for integral to unity
      double norm = dist[gin][i_gout][order - 1];
      if (norm > 0.) {
        for (int imu = 0; imu < order; imu++) {
          fmu[gin][i_gout][imu] /= norm;
          dist[gin][i_gout][imu] /= norm;
        }
      }
    }
  }
}


double ScattDataHistogram::calc_f(int gin, int gout, double mu)
{
  // TODO: gout >= or gout >?
  double f;
  if ((gout < gmin[gin]) || (gout >= gmax[gin])) {
    f = 0.;
  } else {
    // Find mu bin
    int i_gout = gout - gmin[gin]; //TODO: + 1?
    int imu;
    if (mu == 1.) {
      // use size -2 to have the index one before the end
      imu = this->mu.size() - 2;
    } else {
      imu = std::floor((mu + 1.) / dmu + 1.);
    }

    f = fmu[gin][i_gout][imu];
  }
  return f;
}


void ScattDataHistogram::sample(int gin, int& gout, double& mu, double& wgt)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout);

  // Determine the outgoing cosine bin
  double xi = prn();

  int imu;
  if (xi < dist[gin][i_gout][0]) {
    imu = 0;
  } else {
    // TODO lower_bound?  + 1?
    imu = std::upper_bound(dist[gin][i_gout].begin(),
                           dist[gin][i_gout].end(), xi) -
         dist[gin][i_gout].begin();
  }

  // Randomly select mu within the imu bin
  mu = prn() * dmu + this->mu[imu];

  if (mu < -1.) {
    mu = -1.;
  } else if (mu > 1.) {
    mu = 1.;
  }

  // Update the weight to reflect neutron multiplicity
  wgt *= mult[gin][i_gout];
}


double_3dvec ScattDataHistogram::get_matrix(int max_order)
{
  // Get the sizes and initialize the data to 0
  int groups = energy.size();
  // We ignore the requested order for Histogram and Tabular representations
  int order_dim = get_order();
  double_3dvec matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(order_dim, 0.)));

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix[gin][gout][l] = scattxs[gin] * energy[gin][i_gout] *
             fmu[gin][i_gout][l];
      }
    }
  }
  return matrix;
}


void ScattDataHistogram::combine(std::vector<ScattData*>& those_scatts,
                                double_1dvec& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  int max_order;
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataHistogram* that = dynamic_cast<ScattDataHistogram*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    if (i == 0) {
      max_order = that->get_order();
    } else if (max_order != that->get_order()) {
      fatal_error("Cannot combine the ScattData objects!");
    }
  }

  // Get the groups as a shorthand
  int groups = dynamic_cast<ScattDataHistogram*>(those_scatts[0])->energy.size();

  // Now allocate and zero our storage spaces
  double_3dvec this_matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(max_order, 0.)));
  double_2dvec mult_numer(groups, double_1dvec(groups, 0.));
  double_2dvec mult_denom(groups, double_1dvec(groups, 0.));

  // Build the dense scattering and multiplicity matrices
  // Get the multiplicity_matrix
  // To combine from nuclidic data we need to use the final relationship
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // Developed as follows:
  // mult_{gg'} = nuScatt{g,g'} / Scatt{g,g'}
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) / sum(N_i*scatt_{i,g,g'})
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // nuscatt_{i,g,g'} can be reconstructed from the energy and scattxs member
  // variables
  for (int i = 0; i < those_scatts.size(); i++) {
    ScattDataHistogram* that = dynamic_cast<ScattDataHistogram*>(those_scatts[i]);

    // Build the dense matrix for that object
    double_3dvec that_matrix = that->get_matrix(max_order);

    // Now add that to this for the scattering and multiplicity
    for (int gin = 0; gin < groups; gin++) {
      // Only spend time adding that's gmin to gmax data since the rest will
      // be zeros
      int i_gout = 0;
      for (int gout = that->gmin[gin]; gout <= that->gmax[gin]; gout++) {
        // Do the scattering matrix
        for (int l = 0; l < max_order; l++) {
          this_matrix[gin][gout][l] += scalars[i] * that_matrix[gin][gout][l];
        }

        // Incorporate that's contribution to the multiplicity matrix data
        double nuscatt = that->scattxs[gin] * that->energy[gin][i_gout];
        mult_numer[gin][gout] += scalars[i] * nuscatt;
        if (that->mult[gin][i_gout] > 0.) {
          mult_denom[gin][gout] += scalars[i] * nuscatt / that->mult[gin][i_gout];
        } else {
          mult_denom[gin][gout] += scalars[i];
        }
        i_gout++;
      }
    }
  }

  // Combine mult_numer and mult_denom into the combined multiplicity matrix
  double_2dvec this_mult(groups, double_1dvec(groups, 1.));
  for (int gin = 0; gin < groups; gin++) {
    for (int gout = 0; gout < groups; gout++) {
      if (mult_denom[gin][gout] > 0.) {
        this_mult[gin][gout] = mult_numer[gin][gout] / mult_denom[gin][gout];
      }
    }
  }
  mult_numer.clear();
  mult_denom.clear();

  // We have the data, now we need to convert to a jagged array and then use
  // the initialize function to store it on the object.
  int_1dvec in_gmin(groups);
  int_1dvec in_gmax(groups);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);
  for (int gin = 0; gin < groups; gin++) {
    // Find the minimum and maximum group boundaries
    int gmin_;
    for (gmin_ = 0; gmin_ < groups; gmin_++) {
      bool non_zero = std::all_of(this_matrix[gin][gmin_].begin(),
                                  this_matrix[gin][gmin_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }
    int gmax_;
    for (gmax_ = groups - 1; gmax_ >= 0; gmax_--) {
      bool non_zero = std::all_of(this_matrix[gin][gmax_].begin(),
                                  this_matrix[gin][gmax_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }

    // treat the case of all values being 0
    if (gmin_ > gmax_) {
      gmin_ = gin;
      gmax_ = gin;
    }

    // Store the group bounds
    in_gmin[gin] = gmin_;
    in_gmax[gin] = gmax_;

    // Store the data in the compressed format
    sparse_scatter[gin].resize(gmax_ - gmin_ + 1);
    sparse_mult[gin].resize(gmax_ - gmin_ + 1);
    int i_gout = 0;
    for (int gout = gmin_; gout <= gmax_; gout++) {
      sparse_scatter[gin][i_gout] = this_matrix[gin][gout];
      sparse_mult[gin][i_gout] = this_mult[gin][gout];
    }
  }

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}

//==============================================================================
// ScattDataTabular methods
//==============================================================================

void ScattDataTabular::init(int_1dvec& in_gmin, int_1dvec& in_gmax,
                            double_2dvec& in_mult, double_3dvec& coeffs)
{
  int groups = coeffs.size();
  int order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Build the angular distribution mu values
  mu = double_1dvec(order);
  dmu = 2. / (order - 1);
  mu[0] = -1.;
  for (int imu = 1; imu < order - 1; imu++) {
    mu[imu] = -1. + imu * dmu;
  }
  mu[order - 1] = 1.;

  // Get the scattering cross section value by integrating the distribution
  // over all mu points and then combining over all outgoing groups
  scattxs.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    scattxs[gin] = 0.;
    for (int i_gout = 0; i_gout < matrix[gin].size(); i_gout++) {
      for (int imu = 1; imu < order; imu++) {
        scattxs[gin] += 0.5 * dmu * (matrix[gin][i_gout][imu - 1] +
                                     matrix[gin][i_gout][imu]);
      }
    }
  }

  // Build the energy transfer matrix from data in the variable matrix
  double_2dvec in_energy;
  in_energy.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = in_gmax[gin] - in_gmin[gin] + 1;
    in_energy[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      double norm = 0.;
      for (int imu = 1; imu < order; imu++) {
        norm += 0.5 * dmu * (matrix[gin][i_gout][imu - 1] +
                             matrix[gin][i_gout][imu]);
      }
      in_energy[gin][i_gout] = norm;
    }
  }

  // Initialize the base class attributes
  ScattData::generic_init(order, in_gmin, in_gmax, in_energy, in_mult);

  // Calculate f(mu) and integrate it so we can avoid rejection sampling
  fmu.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = gmax[gin] - gmin[gin] + 1;
    fmu[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      fmu[gin][i_gout].resize(order);
      // The variable matrix contains f(mu); so directly assign it
      fmu[gin][i_gout] = matrix[gin][i_gout];

      // Ensure positivity
      for (auto& val : fmu[gin][i_gout]) {
        if (val < 0.) val = 0.;
      }

      // Now re-normalize for numerical integration issues and to take care of
      // the above negative fix-up.  Also accrue the CDF
      double norm = 0.;
      for (int imu = 1; imu < order; imu++) {
        norm += 0.5 * dmu * (fmu[gin][i_gout][imu - 1] +
                             fmu[gin][i_gout][imu]);
        // incorporate to the CDF
        dist[gin][i_gout][imu] = norm;
      }

      // now do the normalization
      if (norm > 0.) {
        for (int imu = 0; imu < order; imu++) {
          fmu[gin][i_gout][imu] /= norm;
          dist[gin][i_gout][imu] /= norm;
        }
      }
    }
  }
}


double ScattDataTabular::calc_f(int gin, int gout, double mu)
{
  // TODO: gout >= or gout >?
  double f;
  if ((gout < gmin[gin]) || (gout >= gmax[gin])) {
    f = 0.;
  } else {
    // Find mu bin
    int i_gout = gout - gmin[gin]; //TODO: + 1?
    int imu;
    if (mu == 1.) {
      // use size -2 to have the index one before the end
      imu = this->mu.size() - 2;
    } else {
      imu = std::floor((mu + 1.) / dmu + 1.);
    }

    double r = (mu - this->mu[imu]) / (this->mu[imu + 1] - this->mu[imu]);
    f = (1. - r) * fmu[gin][i_gout][imu] + r * fmu[gin][i_gout][imu + 1];
  }
  return f;
}


void ScattDataTabular::sample(int gin, int& gout, double& mu, double& wgt)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout);

  // Determine the outgoing cosine bin
  int NP = this->mu.size();
  double xi = prn();

  double c_k = dist[gin][i_gout][0];
  int k;
  for (k = 0; k < NP - 1; k++) {
    double c_k1 = dist[gin][i_gout][k + 1];
    if (xi < c_k1) break;
    c_k = c_k1;
  }

  // Check to make sure k is <= NP - 1
  k = std::min(k, NP - 2);

  // Find the pdf values we want
  double p0 = fmu[gin][i_gout][k];
  double mu0 = this -> mu[k];
  double p1 = fmu[gin][i_gout][k + 1];
  double mu1 = this -> mu[k + 1];

  if (p0 == p1) {
    mu = mu0 + (xi - c_k) / p0;
  } else {
    double frac = (p1 - p0) / (mu1 - mu0);
    mu = mu0 + (std::sqrt(std::max(0., p0 * p0 + 2. * frac * (xi - c_k)))
                - p0) / frac;
  }

  if (mu < -1.) {
    mu = -1.;
  } else if (mu > 1.) {
    mu = 1.;
  }

  // Update the weight to reflect neutron multiplicity
  wgt *= mult[gin][i_gout];
}


double_3dvec ScattDataTabular::get_matrix(int max_order)
{
  // Get the sizes and initialize the data to 0
  int groups = energy.size();
  // We ignore the requested order for Histogram and Tabular representations
  int order_dim = get_order();
  double_3dvec matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(order_dim, 0.)));

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix[gin][gout][l] = scattxs[gin] * energy[gin][i_gout] *
             fmu[gin][i_gout][l];
      }
    }
  }
  return matrix;
}

void ScattDataTabular::combine(std::vector<ScattData*>& those_scatts,
                               double_1dvec& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  int max_order;
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataTabular* that = dynamic_cast<ScattDataTabular*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    if (i == 0) {
      max_order = that->get_order();
    } else if (max_order != that->get_order()) {
      fatal_error("Cannot combine the ScattData objects!");
    }
  }

  // Get the groups as a shorthand
  int groups = dynamic_cast<ScattDataTabular*>(those_scatts[0])->energy.size();

  // Now allocate and zero our storage spaces
  double_3dvec this_matrix = double_3dvec(groups, double_2dvec(groups,
       double_1dvec(max_order, 0.)));
  double_2dvec mult_numer(groups, double_1dvec(groups, 0.));
  double_2dvec mult_denom(groups, double_1dvec(groups, 0.));

  // Build the dense scattering and multiplicity matrices
  // Get the multiplicity_matrix
  // To combine from nuclidic data we need to use the final relationship
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // Developed as follows:
  // mult_{gg'} = nuScatt{g,g'} / Scatt{g,g'}
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) / sum(N_i*scatt_{i,g,g'})
  // mult_{gg'} = sum_i(N_i*nuscatt_{i,g,g'}) /
  //              sum_i(N_i*(nuscatt_{i,g,g'} / mult_{i,g,g'}))
  // nuscatt_{i,g,g'} can be reconstructed from the energy and scattxs member
  // variables
  for (int i = 0; i < those_scatts.size(); i++) {
    ScattDataTabular* that = dynamic_cast<ScattDataTabular*>(those_scatts[i]);

    // Build the dense matrix for that object
    double_3dvec that_matrix = that->get_matrix(max_order);

    // Now add that to this for the scattering and multiplicity
    for (int gin = 0; gin < groups; gin++) {
      // Only spend time adding that's gmin to gmax data since the rest will
      // be zeros
      int i_gout = 0;
      for (int gout = that->gmin[gin]; gout <= that->gmax[gin]; gout++) {
        // Do the scattering matrix
        for (int l = 0; l < max_order; l++) {
          this_matrix[gin][gout][l] += scalars[i] * that_matrix[gin][gout][l];
        }

        // Incorporate that's contribution to the multiplicity matrix data
        double nuscatt = that->scattxs[gin] * that->energy[gin][i_gout];
        mult_numer[gin][gout] += scalars[i] * nuscatt;
        if (that->mult[gin][i_gout] > 0.) {
          mult_denom[gin][gout] += scalars[i] * nuscatt / that->mult[gin][i_gout];
        } else {
          mult_denom[gin][gout] += scalars[i];
        }
        i_gout++;
      }
    }
  }

  // Combine mult_numer and mult_denom into the combined multiplicity matrix
  double_2dvec this_mult(groups, double_1dvec(groups, 1.));
  for (int gin = 0; gin < groups; gin++) {
    for (int gout = 0; gout < groups; gout++) {
      if (mult_denom[gin][gout] > 0.) {
        this_mult[gin][gout] = mult_numer[gin][gout] / mult_denom[gin][gout];
      }
    }
  }
  mult_numer.clear();
  mult_denom.clear();

  // We have the data, now we need to convert to a jagged array and then use
  // the initialize function to store it on the object.
  int_1dvec in_gmin(groups);
  int_1dvec in_gmax(groups);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);
  for (int gin = 0; gin < groups; gin++) {
    // Find the minimum and maximum group boundaries
    int gmin_;
    for (gmin_ = 0; gmin_ < groups; gmin_++) {
      bool non_zero = std::all_of(this_matrix[gin][gmin_].begin(),
                                  this_matrix[gin][gmin_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }
    int gmax_;
    for (gmax_ = groups - 1; gmax_ >= 0; gmax_--) {
      bool non_zero = std::all_of(this_matrix[gin][gmax_].begin(),
                                  this_matrix[gin][gmax_].end(),
                                  [](double val){return val != 0.;});
      if (non_zero) break;
    }

    // treat the case of all values being 0
    if (gmin_ > gmax_) {
      gmin_ = gin;
      gmax_ = gin;
    }

    // Store the group bounds
    in_gmin[gin] = gmin_;
    in_gmax[gin] = gmax_;

    // Store the data in the compressed format
    sparse_scatter[gin].resize(gmax_ - gmin_ + 1);
    sparse_mult[gin].resize(gmax_ - gmin_ + 1);
    int i_gout = 0;
    for (int gout = gmin_; gout <= gmax_; gout++) {
      sparse_scatter[gin][i_gout] = this_matrix[gin][gout];
      sparse_mult[gin][i_gout] = this_mult[gin][gout];
      i_gout++;
    }
  }

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}


void convert_legendre_to_tabular(ScattDataLegendre& leg,
                                 ScattDataTabular& tab, int n_mu)
{
  tab.generic_init(n_mu, leg.gmin, leg.gmax, leg.energy, leg.mult);
  tab.scattxs = leg.scattxs;

  // Build mu and dmu
  tab.mu = double_1dvec(n_mu);
  tab.dmu = 2. / (n_mu - 1);
  tab.mu[0] = -1.;
  for (int imu = 1; imu < n_mu - 1; imu++) {
    tab.mu[imu] = -1. + (imu - 1) * tab.dmu;
  }
  tab.mu[n_mu - 1] = 1.;

  // Calculate f(mu) and integrate it so we can avoid rejection sampling
  int groups = tab.energy.size();
  tab.fmu.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = tab.gmax[gin] - tab.gmin[gin] + 1;
    tab.fmu[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      tab.fmu[gin][i_gout].resize(n_mu);
      for (int imu = 0; imu < n_mu; imu++) {
        tab.fmu[gin][i_gout][imu] =
             evaluate_legendre_c(leg.dist[gin][i_gout].size() - 1,
                                 leg.dist[gin][i_gout].data(), tab.mu[imu]);
      }

      // Ensure positivity
      for (auto& val : tab.fmu[gin][i_gout]) {
        if (val < 0.) val = 0.;
      }

      // Now re-normalize for numerical integration issues and to take care of
      // the above negative fix-up.  Also accrue the CDF
      double norm = 0.;
      for (int imu = 1; imu < n_mu; imu++) {
        norm += 0.5 * tab.dmu * (tab.fmu[gin][i_gout][imu - 1] +
                                 tab.fmu[gin][i_gout][imu]);
        // incorporate to the CDF
        tab.dist[gin][i_gout][imu] = norm;
      }

      // now do the normalization
      if (norm > 0.) {
        for (int imu = 0; imu < n_mu; imu++) {
          tab.fmu[gin][i_gout][imu] /= norm;
          tab.dist[gin][i_gout][imu] /= norm;
        }
      }
    }
  }
}

} // namespace openmc
