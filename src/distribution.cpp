#include "distribution.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "error.h"
#include "math_functions.h"
#include "random_lcg.h"

namespace openmc {

Discrete::Discrete(const double* x, const double* p, int n)
  : x_{x, x+n}, p_{p, p+n}
{
  // Renormalize density function so that it sums to unity
  double norm = std::accumulate(p_.begin(), p_.end(), 0.0);
  for (auto& p_i : p_)
    p_i /= norm;
}


double Discrete::sample()
{
  int n = x_.size();
  if (n > 1) {
    double xi = prn();
    double c = 0.0;
    for (int i = 0; i < n; ++i) {
      c += p_[i];
      if (xi < c) return x_[i];
    }
    // throw exception?
  } else {
    return x_[0];
  }
}


double Uniform::sample()
{
  return a_ + prn()*(b_ - a_);
}


double Maxwell::sample()
{
  return maxwell_spectrum_c(theta_);
}


double Watt::sample()
{
  return watt_spectrum_c(a_, b_);
}


Tabular::Tabular(const double* x, const double* p, int n, Interpolation interp)
  : x_{x, x+n}, p_{p, p+n}, interp_{interp}, c_(n, 0.0)
{
  // Check interpolation parameter
  if (interp_ != Interpolation::histogram &&
      interp_ != Interpolation::lin_lin) {
    openmc::fatal_error("Only histogram and linear-linear interpolation "
                        "for tabular distribution is supported.");
  }

  // Calculate cumulative distribution function
  for (int i = 1; i < n; ++i) {
    if (interp_ == Interpolation::histogram) {
      c_[i] = c_[i-1] + p_[i-1]*(x_[i] - x_[i-1]);
    } else if (interp_ == Interpolation::lin_lin) {
      c_[i] = c_[i-1] + 0.5*(p_[i-1] + p_[i]) * (x_[i] - x_[i-1]);
    }
  }

  // Normalize density and distribution functions
  for (int i = 0; i < n; ++i) {
    p_[i] = p_[i]/c_[n-1];
    c_[i] = c_[i]/c_[n-1];
  }
}


double Tabular::sample()
{
  // Sample value of CDF
  double c = prn();

  // Find first CDF bin which is above the sampled value
  double c_i = c_[0];
  int i;
  size_t n = c_.size();
  for (i = 0; i < n - 1; ++i) {
    if (c <= c_[i+1]) break;
    c_i = c_[i+1];
  }

  // Determine bounding PDF values
  double x_i = x_[i];
  double p_i = p_[i];

  if (interp_ == Interpolation::histogram) {
    // Histogram interpolation
    if (p_i > 0.0) {
      return x_i + (c - c_i)/p_i;
    } else {
      return x_i;
    }
  } else {
    // Linear-linear interpolation
    double x_i1 = x_[i + 1];
    double p_i1 = p_[i + 1];

    double m = (p_i1 - p_i)/(x_i1 - x_i);
    if (m == 0.0) {
      return x_i + (c - c_i)/p_i;
    } else {
      return x_i + (std::sqrt(std::max(0.0, p_i*p_i + 2*m*(c - c_i))) - p_i)/m;
    }
  }
}


double Equiprobable::sample()
{
  size_t n = x_.size();

  double r = prn();
  int i = std::floor((n - 1)*r);

  double xl = x_[i];
  double xr = x_[i+i];
  return xl + ((n - 1)*r - i) * (xr - xl);
}

} // namespace openmc
