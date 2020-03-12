#include "openmc/scattdata.h"

#include <algorithm>
#include <numeric>
#include <cmath>

#include "xtensor/xbuilder.hpp"

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// ScattData base-class methods
//==============================================================================

void
ScattData::base_init(int order, const xt::xtensor<int, 1>& in_gmin,
     const xt::xtensor<int, 1>& in_gmax, const double_2dvec& in_energy,
     const double_2dvec& in_mult)
{
  size_t groups = in_energy.size();

  gmin = in_gmin;
  gmax = in_gmax;
  energy.resize(groups);
  mult.resize(groups);
  dist.resize(groups);

  for (int gin = 0; gin < groups; gin++) {
    // Store the inputted data
    energy[gin] = in_energy[gin];
    mult[gin] = in_mult[gin];

    // Make sure the energy is normalized
    double norm = std::accumulate(energy[gin].begin(), energy[gin].end(), 0.);

    if (norm != 0.) {
      for (auto& n : energy[gin]) n /= norm;
    }

    // Initialize the distribution data
    dist[gin].resize(in_gmax[gin] - in_gmin[gin] + 1);
    for (auto& v : dist[gin]) {
      v.resize(order);
    }
  }
}

//==============================================================================

void
ScattData::base_combine(size_t max_order,
     const std::vector<ScattData*>& those_scatts, const std::vector<double>& scalars,
     xt::xtensor<int, 1>& in_gmin, xt::xtensor<int, 1>& in_gmax, double_2dvec& sparse_mult,
     double_3dvec& sparse_scatter)
{
  size_t groups = those_scatts[0] -> energy.size();

  // Now allocate and zero our storage spaces
  xt::xtensor<double, 3> this_matrix({groups, groups, max_order}, 0.);
  xt::xtensor<double, 2> mult_numer({groups, groups}, 0.);
  xt::xtensor<double, 2> mult_denom({groups, groups}, 0.);
  // TODO: Need to review this:
  if (this->scattxs.size() > 0) {
    this_matrix = this->get_matrix(max_order);
  }
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
    ScattData* that = those_scatts[i];

    // Build the dense matrix for that object
    xt::xtensor<double, 3> that_matrix = that->get_matrix(max_order);

    // Now add that to this for the scattering and multiplicity
    for (int gin = 0; gin < groups; gin++) {
      // Only spend time adding that's gmin to gmax data since the rest will
      // be zeros
      int i_gout = 0;
      for (int gout = that->gmin(gin); gout <= that->gmax(gin); gout++) {
        // Do the scattering matrix
        for (int l = 0; l < max_order; l++) {
          this_matrix(gin, gout, l) += scalars[i] * that_matrix(gin, gout, l);
        }

        // Incorporate that's contribution to the multiplicity matrix data
        double nuscatt = that->scattxs(gin) * that->energy[gin][i_gout];
        mult_numer(gin, gout) += scalars[i] * nuscatt;
        if (that->mult[gin][i_gout] > 0.) {
          mult_denom(gin, gout) += scalars[i] * nuscatt / that->mult[gin][i_gout];
        } else {
          mult_denom(gin, gout) += scalars[i];
        }
        i_gout++;
      }
    }
  }

  // Combine mult_numer and mult_denom into the combined multiplicity matrix
  xt::xtensor<double, 2> this_mult({groups, groups}, 1.);
  // TODO: Need to check this too
  for (int gin = 0; gin < groups; gin++) {
    for (int gout = 0; gout < groups; gout++) {
      if (std::abs(mult_denom(gin, gout)) > 0.0) {
	   this_mult(gin, gout) = mult_numer(gin, gout) / mult_denom(gin, gout);
      } else {
	   if (mult_numer(gin, gout) == 0.0) {
	     this_mult(gin, gout) = 1.0;
	   }
      }
    }
  }

  // We have the data, now we need to convert to a jagged array and then use
  // the initialize function to store it on the object.
  for (int gin = 0; gin < groups; gin++) {
    // Find the minimum and maximum group boundaries
    int gmin_;
    for (gmin_ = 0; gmin_ < groups; gmin_++) {
      bool non_zero = false;
      for (int l = 0; l < this_matrix.shape()[2]; l++) {
        if (this_matrix(gin, gmin_, l) != 0.) {
          non_zero = true;
          break;
        }
      }
      if (non_zero) break;
    }
    int gmax_;
    for (gmax_ = groups - 1; gmax_ >= 0; gmax_--) {
      bool non_zero = false;
      for (int l = 0; l < this_matrix.shape()[2]; l++) {
        if (this_matrix(gin, gmax_, l) != 0.) {
          non_zero = true;
          break;
        }
      }
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
      sparse_scatter[gin][i_gout].resize(this_matrix.shape()[2]);
      for (int l = 0; l < this_matrix.shape()[2]; l++) {
        sparse_scatter[gin][i_gout][l] = this_matrix(gin, gout, l);
      }
      sparse_mult[gin][i_gout] = this_mult(gin, gout);
      i_gout++;
    }
  }
}

//==============================================================================

void
ScattData::sample_energy(int gin, int& gout, int& i_gout, uint64_t* seed)
{
  // Sample the outgoing group
  double xi = prn(seed);
  double prob = 0.;
  i_gout = 0;
  for (gout = gmin[gin]; gout < gmax[gin]; ++gout) {
    prob += energy[gin][i_gout];
    if (xi < prob) break;
    ++i_gout;
  }
}

//==============================================================================

double
ScattData::get_xs(MgxsType xstype, int gin, const int* gout, const double* mu)
{
  // Set the outgoing group offset index as needed
  int i_gout = 0;
  if (gout != nullptr) {
    // short circuit the function if gout is from a zero portion of the
    // scattering matrix
    if ((*gout < gmin[gin]) || (*gout > gmax[gin])) { // > gmax?
      return 0.;
    }
    i_gout = *gout - gmin[gin];
  }

  double val = scattxs[gin];
  switch(xstype) {
  case MgxsType::SCATTER:
    if (gout != nullptr) val *= energy[gin][i_gout];
    break;
  case MgxsType::SCATTER_MULT:
    if (gout != nullptr) {
      val *= energy[gin][i_gout] / mult[gin][i_gout];
    } else {
      val /= std::inner_product(mult[gin].begin(), mult[gin].end(),
                                energy[gin].begin(), 0.0);
    }
    break;
  case MgxsType::SCATTER_FMU_MULT:
    if ((gout != nullptr) && (mu != nullptr)) {
      val *= energy[gin][i_gout] * calc_f(gin, *gout, *mu);
    } else {
      // This is not an expected path (asking for f_mu without asking for a
      // group or mu is not useful
      fatal_error("Invalid call to get_xs");
    }
    break;
  case MgxsType::SCATTER_FMU:
    if ((gout != nullptr) && (mu != nullptr)) {
      val *= energy[gin][i_gout] * calc_f(gin, *gout, *mu) / mult[gin][i_gout];
    } else {
      // This is not an expected path (asking for f_mu without asking for a
      // group or mu is not useful
      fatal_error("Invalid call to get_xs");
    }
    break;
  default:
    break;
  }
  return val;
}

//==============================================================================
// ScattDataLegendre methods
//==============================================================================

void
ScattDataLegendre::init(const xt::xtensor<int, 1>& in_gmin,
     const xt::xtensor<int, 1>& in_gmax, const double_2dvec& in_mult,
     const double_3dvec& coeffs)
{
  size_t groups = coeffs.size();
  size_t order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Get the scattering cross section value by summing the un-normalized P0
  // coefficient in the variable matrix over all outgoing groups.
  scattxs = xt::zeros<double>({groups});
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = in_gmax[gin] - in_gmin[gin] + 1;
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      scattxs[gin] += matrix[gin][i_gout][0];
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
  ScattData::base_init(order, in_gmin, in_gmax, in_energy, in_mult);

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

//==============================================================================

void
ScattDataLegendre::update_max_val()
{
  size_t groups = max_val.size();
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
        double f = evaluate_legendre(dist[gin][i_gout].size() - 1,
                                     dist[gin][i_gout].data(), mu);

        // if this is a new maximum, store it
        if (f > max_val[gin][i_gout]) max_val[gin][i_gout] = f;
      } // end imu loop

      // Since we may not have caught the true max, add 10% margin
      max_val[gin][i_gout] *= 1.1;
    }
  }
}

//==============================================================================

double
ScattDataLegendre::calc_f(int gin, int gout, double mu)
{
  double f;
  if ((gout < gmin[gin]) || (gout > gmax[gin])) {
    f = 0.;
  } else {
    int i_gout = gout - gmin[gin];
    f = evaluate_legendre(dist[gin][i_gout].size() - 1,
                          dist[gin][i_gout].data(), mu);
  }
  return f;
}

//==============================================================================

void
ScattDataLegendre::sample(int gin, int& gout, double& mu, double& wgt,
                          uint64_t* seed)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout, seed);

  // Now we can sample mu using the scattering kernel using rejection
  // sampling from a rectangular bounding box
  double M = max_val[gin][i_gout];
  int samples;
  for (samples = 0; samples < MAX_SAMPLE; ++samples) {
    mu = 2. * prn(seed) - 1.;
    double f = calc_f(gin, gout, mu);
    if (f > 0.) {
      double u = prn(seed) * M;
      if (u <= f) break;
    }
  }
  if (samples == MAX_SAMPLE) {
    fatal_error("Maximum number of Legendre expansion samples reached!");
  }

  // Update the weight to reflect neutron multiplicity
  wgt *= mult[gin][i_gout];
}

//==============================================================================

void
ScattDataLegendre::combine(const std::vector<ScattData*>& those_scatts,
                           const std::vector<double>& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  size_t max_order = 0;
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataLegendre* that = dynamic_cast<ScattDataLegendre*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    size_t that_order = that->get_order();
    if (that_order > max_order) max_order = that_order;
  }
  max_order++;  // Add one since this is a Legendre

  size_t groups = those_scatts[0] -> energy.size();

  xt::xtensor<int, 1> in_gmin({groups}, 0);
  xt::xtensor<int, 1> in_gmax({groups}, 0);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);

  // The rest of the steps do not depend on the type of angular representation
  // so we use a base class method to sum up xs and create new energy and mult
  // matrices
  ScattData::base_combine(max_order, those_scatts, scalars, in_gmin, in_gmax,
                          sparse_mult, sparse_scatter);

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}

//==============================================================================

xt::xtensor<double, 3>
ScattDataLegendre::get_matrix(size_t max_order)
{
  // Get the sizes and initialize the data to 0
  size_t groups = energy.size();
  size_t order_dim = max_order + 1;
  xt::xtensor<double, 3> matrix({groups, groups, order_dim}, 0.);

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix(gin, gout, l) = scattxs[gin] * energy[gin][i_gout] *
             dist[gin][i_gout][l];
      }
    }
  }
  return matrix;
}

//==============================================================================
// ScattDataHistogram methods
//==============================================================================

void
ScattDataHistogram::init(const xt::xtensor<int, 1>& in_gmin,
     const xt::xtensor<int, 1>& in_gmax, const double_2dvec& in_mult,
     const double_3dvec& coeffs)
{
  size_t groups = coeffs.size();
  size_t order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Get the scattering cross section value by summing the distribution
  // over all the histogram bins in angle and outgoing energy groups
  scattxs = xt::zeros<double>({groups});
  for (int gin = 0; gin < groups; gin++) {
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
  ScattData::base_init(order, in_gmin, in_gmax, in_energy, in_mult);

  // Build the angular distribution mu values
  mu = xt::linspace(-1., 1., order + 1);
  dmu = 2. / order;

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

//==============================================================================

double
ScattDataHistogram::calc_f(int gin, int gout, double mu)
{
  double f;
  if ((gout < gmin[gin]) || (gout > gmax[gin])) {
    f = 0.;
  } else {
    // Find mu bin
    int i_gout = gout - gmin[gin];
    int imu;
    if (mu == 1.) {
      // use size -2 to have the index one before the end
      imu = this->mu.shape()[0] - 2;
    } else {
      imu = std::floor((mu + 1.) / dmu + 1.) - 1;
    }

    f = fmu[gin][i_gout][imu];
  }
  return f;
}

//==============================================================================

void
ScattDataHistogram::sample(int gin, int& gout, double& mu, double& wgt,
                           uint64_t* seed)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout, seed);

  // Determine the outgoing cosine bin
  double xi = prn(seed);

  int imu;
  if (xi < dist[gin][i_gout][0]) {
    imu = 0;
  } else {
    imu = std::upper_bound(dist[gin][i_gout].begin(),
                           dist[gin][i_gout].end(), xi) -
         dist[gin][i_gout].begin();
  }

  // Randomly select mu within the imu bin
  mu = prn(seed) * dmu + this->mu[imu];

  if (mu < -1.) {
    mu = -1.;
  } else if (mu > 1.) {
    mu = 1.;
  }

  // Update the weight to reflect neutron multiplicity
  wgt *= mult[gin][i_gout];
}

//==============================================================================

xt::xtensor<double, 3>
ScattDataHistogram::get_matrix(size_t max_order)
{
  // Get the sizes and initialize the data to 0
  size_t groups = energy.size();
  // We ignore the requested order for Histogram and Tabular representations
  size_t order_dim = get_order();
  xt::xtensor<double, 3> matrix({groups, groups, order_dim}, 0);

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix(gin, gout, l) = scattxs[gin] * energy[gin][i_gout] *
             fmu[gin][i_gout][l];
      }
    }
  }
  return matrix;
}

//==============================================================================

void
ScattDataHistogram::combine(const std::vector<ScattData*>& those_scatts,
                            const std::vector<double>& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  size_t max_order = those_scatts[0]->get_order();
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataHistogram* that = dynamic_cast<ScattDataHistogram*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    if (max_order != that->get_order()) {
      fatal_error("Cannot combine the ScattData objects!");
    }
  }

  size_t groups = those_scatts[0] -> energy.size();

  xt::xtensor<int, 1> in_gmin({groups}, 0);
  xt::xtensor<int, 1> in_gmax({groups}, 0);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);

  // The rest of the steps do not depend on the type of angular representation
  // so we use a base class method to sum up xs and create new energy and mult
  // matrices
  ScattData::base_combine(max_order, those_scatts, scalars, in_gmin, in_gmax,
                          sparse_mult, sparse_scatter);

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}

//==============================================================================
// ScattDataTabular methods
//==============================================================================

void
ScattDataTabular::init(const xt::xtensor<int, 1>& in_gmin,
     const xt::xtensor<int, 1>& in_gmax, const double_2dvec& in_mult,
     const double_3dvec& coeffs)
{
  size_t groups = coeffs.size();
  size_t order = coeffs[0][0].size();

  // make a copy of coeffs that we can use to both extract data and normalize
  double_3dvec matrix = coeffs;

  // Build the angular distribution mu values
  mu = xt::linspace(-1., 1., order);
  dmu = 2. / (order - 1);

  // Get the scattering cross section value by integrating the distribution
  // over all mu points and then combining over all outgoing groups
  scattxs = xt::zeros<double>({groups});
  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < matrix[gin].size(); i_gout++) {
      for (int imu = 1; imu < order; imu++) {
        scattxs[gin] += 0.5 * dmu * (matrix[gin][i_gout][imu - 1] +
                                     matrix[gin][i_gout][imu]);
      }
    }
  }

  // Build the energy transfer matrix from data in the variable matrix
  double_2dvec in_energy(groups);
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
  ScattData::base_init(order, in_gmin, in_gmax, in_energy, in_mult);

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

//==============================================================================

double
ScattDataTabular::calc_f(int gin, int gout, double mu)
{
  double f;
  if ((gout < gmin[gin]) || (gout > gmax[gin])) {
    f = 0.;
  } else {
    // Find mu bin
    int i_gout = gout - gmin[gin];
    int imu;
    if (mu == 1.) {
      // use size -2 to have the index one before the end
      imu = this->mu.shape()[0] - 2;
    } else {
      imu = std::floor((mu + 1.) / dmu + 1.) - 1;
    }

    double r = (mu - this->mu[imu]) / (this->mu[imu + 1] - this->mu[imu]);
    f = (1. - r) * fmu[gin][i_gout][imu] + r * fmu[gin][i_gout][imu + 1];
  }
  return f;
}

//==============================================================================

void
ScattDataTabular::sample(int gin, int& gout, double& mu, double& wgt,
                         uint64_t* seed)
{
  // Sample the outgoing energy using the base-class method
  int i_gout;
  sample_energy(gin, gout, i_gout, seed);

  // Determine the outgoing cosine bin
  int NP = this->mu.shape()[0];
  double xi = prn(seed);

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

//==============================================================================

xt::xtensor<double, 3>
ScattDataTabular::get_matrix(size_t max_order)
{
  // Get the sizes and initialize the data to 0
  size_t groups = energy.size();
  // We ignore the requested order for Histogram and Tabular representations
  size_t order_dim = get_order();
  xt::xtensor<double, 3> matrix({groups, groups, order_dim}, 0.);

  for (int gin = 0; gin < groups; gin++) {
    for (int i_gout = 0; i_gout < energy[gin].size(); i_gout++) {
      int gout = i_gout + gmin[gin];
      for (int l = 0; l < order_dim; l++) {
        matrix(gin, gout, l) = scattxs[gin] * energy[gin][i_gout] *
             fmu[gin][i_gout][l];
      }
    }
  }
  return matrix;
}

//==============================================================================

void
ScattDataTabular::combine(const std::vector<ScattData*>& those_scatts,
                          const std::vector<double>& scalars)
{
  // Find the max order in the data set and make sure we can combine the sets
  size_t max_order = those_scatts[0]->get_order();
  for (int i = 0; i < those_scatts.size(); i++) {
    // Lets also make sure these items are combineable
    ScattDataTabular* that = dynamic_cast<ScattDataTabular*>(those_scatts[i]);
    if (!that) {
      fatal_error("Cannot combine the ScattData objects!");
    }
    if (max_order != that->get_order()) {
      fatal_error("Cannot combine the ScattData objects!");
    }
  }

  size_t groups = those_scatts[0] -> energy.size();

  xt::xtensor<int, 1> in_gmin({groups}, 0);
  xt::xtensor<int, 1> in_gmax({groups}, 0);
  double_3dvec sparse_scatter(groups);
  double_2dvec sparse_mult(groups);

  // The rest of the steps do not depend on the type of angular representation
  // so we use a base class method to sum up xs and create new energy and mult
  // matrices
  ScattData::base_combine(max_order, those_scatts, scalars, in_gmin, in_gmax,
                          sparse_mult, sparse_scatter);

  // Got everything we need, store it.
  init(in_gmin, in_gmax, sparse_mult, sparse_scatter);
}

//==============================================================================
// module-level methods
//==============================================================================

void
convert_legendre_to_tabular(ScattDataLegendre& leg, ScattDataTabular& tab)
{
  // See if the user wants us to figure out how many points to use
  int n_mu = settings::legendre_to_tabular_points;
  if (n_mu == C_NONE) {
    // then we will use 2 pts if its P0, or the default if a higher order
    // TODO use an error minimization algorithm that also picks n_mu
    if (leg.get_order() == 0) {
      n_mu = 2;
    } else {
      n_mu = DEFAULT_NMU;
    }
  }

  tab.base_init(n_mu, leg.gmin, leg.gmax, leg.energy, leg.mult);
  tab.scattxs = leg.scattxs;

  // Build mu and dmu
  tab.mu = xt::linspace(-1., 1., n_mu);
  tab.dmu = 2. / (n_mu - 1);

  // Calculate f(mu) and integrate it so we can avoid rejection sampling
  size_t groups = tab.energy.size();
  tab.fmu.resize(groups);
  for (int gin = 0; gin < groups; gin++) {
    int num_groups = tab.gmax[gin] - tab.gmin[gin] + 1;
    tab.fmu[gin].resize(num_groups);
    for (int i_gout = 0; i_gout < num_groups; i_gout++) {
      tab.fmu[gin][i_gout].resize(n_mu);
      for (int imu = 0; imu < n_mu; imu++) {
        tab.fmu[gin][i_gout][imu] =
             evaluate_legendre(leg.dist[gin][i_gout].size() - 1,
                               leg.dist[gin][i_gout].data(), tab.mu[imu]);
      }

      // Ensure positivity
      for (auto& val : tab.fmu[gin][i_gout]) {
        if (val < 0.) val = 0.;
      }

      // Now re-normalize for numerical integration issues and to take care of
      // the above negative fix-up.  Also accrue the CDF
      double norm = 0.;
      tab.dist[gin][i_gout][0] = 0.;
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
