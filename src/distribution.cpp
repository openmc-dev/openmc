#include "openmc/distribution.h"

#include <algorithm> // for copy
#include <cmath>     // for sqrt, floor, max
#include <iterator>  // for back_inserter
#include <numeric>   // for accumulate
#include <stdexcept> // for runtime_error
#include <string>    // for string, stod

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Discrete implementation
//==============================================================================

Discrete::Discrete(pugi::xml_node node)
{
  auto params = get_node_array<xsfloat>(node, "parameters");

  std::size_t n = params.size();
  x_.reserve(n / 2);
  p_.reserve(n / 2);
  std::copy(params.begin(), params.begin() + n/2, std::back_inserter(x_));
  std::copy(params.begin() + n/2, params.end(), std::back_inserter(p_));

  normalize();
}

Discrete::Discrete(const xsfloat* x, const xsfloat* p, int n)
  : x_{x, x+n}, p_{p, p+n}
{
  normalize();
}

xsfloat Discrete::sample(uint64_t* seed) const
{
  int n = x_.size();
  if (n > 1) {
    xsfloat xi = prn(seed);
    xsfloat c = 0.0;
    for (int i = 0; i < n; ++i) {
      c += p_[i];
      if (xi < c) return x_[i];
    }
#ifndef __CUDA_ARCH__
    throw std::runtime_error{"Error when sampling probability mass function."};
#else
    return 0;
#endif
  } else {
    return x_[0];
  }
}

void Discrete::normalize()
{
  // Renormalize density function so that it sums to unity
  xsfloat norm = std::accumulate(p_.begin(), p_.end(), 0.0);
  for (auto& p_i : p_) {
    p_i /= norm;
  }
}

//==============================================================================
// Uniform implementation
//==============================================================================

Uniform::Uniform(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 2) {
    openmc::fatal_error("Uniform distribution must have two "
                        "parameters specified.");
  }

  a_ = params.at(0);
  b_ = params.at(1);
}

xsfloat Uniform::sample(uint64_t* seed) const
{
  return a_ + prn(seed)*(b_ - a_);
}

//==============================================================================
// Maxwell implementation
//==============================================================================

Maxwell::Maxwell(pugi::xml_node node)
{
  theta_ = std::stod(get_node_value(node, "parameters"));
}

xsfloat Maxwell::sample(uint64_t* seed) const
{
  return maxwell_spectrum(theta_, seed);
}

//==============================================================================
// Watt implementation
//==============================================================================

Watt::Watt(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 2)
    openmc::fatal_error("Watt energy distribution must have two "
                        "parameters specified.");

  a_ = params.at(0);
  b_ = params.at(1);
}

xsfloat Watt::sample(uint64_t* seed) const
{
  return watt_spectrum(a_, b_, seed);
}

//==============================================================================
// Normal implementation
//==============================================================================
Normal::Normal(pugi::xml_node node)
{
  auto params = get_node_array<double>(node,"parameters");
  if (params.size() != 2) {
    openmc::fatal_error("Normal energy distribution must have two "
                        "parameters specified.");
  }

  mean_value_ = params.at(0);
  std_dev_ = params.at(1);
}

xsfloat Normal::sample(uint64_t* seed) const
{
  return normal_variate(mean_value_, std_dev_, seed);
}

//==============================================================================
// Muir implementation
//==============================================================================
Muir::Muir(pugi::xml_node node)
{
  auto params = get_node_array<double>(node,"parameters");
  if (params.size() != 3) {
    openmc::fatal_error("Muir energy distribution must have three "
                        "parameters specified.");
  }

  e0_ = params.at(0);
  m_rat_ = params.at(1);
  kt_ = params.at(2);
}

xsfloat Muir::sample(uint64_t* seed) const
{
  return muir_spectrum(e0_, m_rat_, kt_, seed);
}

//==============================================================================
// Tabular implementation
//==============================================================================

Tabular::Tabular(pugi::xml_node node)
{
  if (check_for_node(node, "interpolation")) {
    std::string temp = get_node_value(node, "interpolation");
    if (temp == "histogram") {
      interp_ = Interpolation::histogram;
    } else if (temp == "linear-linear") {
      interp_ = Interpolation::lin_lin;
    } else {
      openmc::fatal_error("Unknown interpolation type for distribution: " + temp);
    }
  } else {
    interp_ = Interpolation::histogram;
  }

  // Read and initialize tabular distribution
  auto params = get_node_array<xsfloat>(node, "parameters");
  std::size_t n = params.size() / 2;
  const xsfloat* x = params.data();
  const xsfloat* p = x + n;
  init(x, p, n);
}

Tabular::Tabular(const xsfloat* x, const xsfloat* p, int n, Interpolation interp, const xsfloat* c)
  : interp_{interp}
{
  init(x, p, n, c);
}

void Tabular::init(const xsfloat* x, const xsfloat* p, std::size_t n, const xsfloat* c)
{
  // Copy x/p arrays into vectors
  x_.reserve(n);
  p_.reserve(n);
  std::copy(x, x + n, std::back_inserter(x_));
  std::copy(p, p + n, std::back_inserter(p_));

  // Check interpolation parameter
  if (interp_ != Interpolation::histogram &&
      interp_ != Interpolation::lin_lin) {
    openmc::fatal_error("Only histogram and linear-linear interpolation "
                        "for tabular distribution is supported.");
  }

  // Calculate cumulative distribution function
  if (c) {
    c_.reserve(n);
    std::copy(c, c + n, std::back_inserter(c_));
  } else {
    c_.resize(n);
    c_[0] = 0.0;
    for (int i = 1; i < n; ++i) {
      if (interp_ == Interpolation::histogram) {
        c_[i] = c_[i-1] + p_[i-1]*(x_[i] - x_[i-1]);
      } else if (interp_ == Interpolation::lin_lin) {
        c_[i] = c_[i-1] + 0.5*(p_[i-1] + p_[i]) * (x_[i] - x_[i-1]);
      }
    }
  }

  // Normalize density and distribution functions
  for (int i = 0; i < n; ++i) {
    p_[i] = p_[i]/c_[n-1];
    c_[i] = c_[i]/c_[n-1];
  }
}

xsfloat Tabular::sample(uint64_t* seed) const
{
  // Sample value of CDF
  xsfloat c = prn(seed);

  // Find first CDF bin which is above the sampled value
  xsfloat c_i = c_[0];
  int i;
  std::size_t n = c_.size();
  for (i = 0; i < n - 1; ++i) {
    if (c <= c_[i+1]) break;
    c_i = c_[i+1];
  }

  // Determine bounding PDF values
  xsfloat const& x_i = x_[i];
  xsfloat const& p_i = p_[i];

  if (interp_ == Interpolation::histogram) {
    // Histogram interpolation
    if (p_i > 0.0) {
      return x_i + (c - c_i)/p_i;
    } else {
      return x_i;
    }
  } else {
    // Linear-linear interpolation
    xsfloat const& x_i1 = x_[i + 1];
    xsfloat const& p_i1 = p_[i + 1];

    xsfloat m = (p_i1 - p_i)/(x_i1 - x_i);
    if (m == 0.0) {
      return x_i + (c - c_i)/p_i;
    } else {
      return x_i + (std::sqrt(std::max(XSZERO, p_i*p_i + 2*m*(c - c_i))) - p_i)/m;
    }
  }
}

//==============================================================================
// Equiprobable implementation
//==============================================================================

xsfloat Equiprobable::sample(uint64_t* seed) const
{
  std::size_t n = x_.size();

  xsfloat r = prn(seed);
  int i = std::floor((n - 1)*r);

  xsfloat const& xl = x_[i];
  xsfloat const& xr = x_[i+i];
  return xl + ((n - 1)*r - i) * (xr - xl);
}

//==============================================================================
// Helper function
//==============================================================================

unique_ptr<Distribution> distribution_from_xml(pugi::xml_node node)
{
  if (!check_for_node(node, "type"))
    openmc::fatal_error("Distribution type must be specified.");

  // Determine type of distribution
  std::string type = get_node_value(node, "type", true, true);

  // Allocate extension of Distribution
  if (type == "uniform") {
    return make_unique<Uniform>(node);
  } else if (type == "maxwell") {
    return make_unique<Maxwell>(node);
  } else if (type == "watt") {
    return make_unique<Watt>(node);
  } else if (type == "normal") {
    return make_unique<Normal>(node);
  } else if (type == "muir") {
    return make_unique<Muir>(node);
  } else if (type == "discrete") {
    return make_unique<Discrete>(node);
  } else if (type == "tabular") {
    return make_unique<Tabular>(node);
  } else {
    openmc::fatal_error("Invalid distribution type: " + type);
    return make_unique<Uniform>(node); // avoid CUDA compiler warning (sad)
  }
}

} // namespace openmc
