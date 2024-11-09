#include "openmc/distribution.h"

#include <algorithm> // for copy
#include <array>
#include <cmath>     // for sqrt, floor, max
#include <iterator>  // for back_inserter
#include <numeric>   // for accumulate
#include <stdexcept> // for runtime_error
#include <string>    // for string, stod

#include <gsl/gsl-lite.hpp>

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// DiscreteIndex implementation
//==============================================================================

DiscreteIndex::DiscreteIndex(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  std::size_t n = params.size() / 2;

  assign({params.data() + n, n});
}

DiscreteIndex::DiscreteIndex(gsl::span<const double> p)
{
  assign(p);
}

void DiscreteIndex::assign(gsl::span<const double> p)
{
  prob_.assign(p.begin(), p.end());

  this->init_alias();
}

void DiscreteIndex::init_alias()
{
  normalize();

  // The initialization and sampling method is based on Vose
  // (DOI: 10.1109/32.92917)
  // Vectors for large and small probabilities based on 1/n
  vector<size_t> large;
  vector<size_t> small;

  size_t n = prob_.size();

  // Set and allocate memory
  alias_.assign(n, 0);

  // Fill large and small vectors based on 1/n
  for (size_t i = 0; i < n; i++) {
    prob_[i] *= n;
    if (prob_[i] > 1.0) {
      large.push_back(i);
    } else {
      small.push_back(i);
    }
  }

  while (!large.empty() && !small.empty()) {
    int j = small.back();
    int k = large.back();

    // Remove last element of small
    small.pop_back();

    // Update probability and alias based on Vose's algorithm
    prob_[k] += prob_[j] - 1.0;
    alias_[j] = k;

    // Move large index to small vector, if it is no longer large
    if (prob_[k] < 1.0) {
      small.push_back(k);
      large.pop_back();
    }
  }
}

size_t DiscreteIndex::sample(uint64_t* seed) const
{
  // Alias sampling of discrete distribution
  size_t n = prob_.size();
  if (n > 1) {
    size_t u = prn(seed) * n;
    if (prn(seed) < prob_[u]) {
      return u;
    } else {
      return alias_[u];
    }
  } else {
    return 0;
  }
}
assign
void DiscreteIndex::normalize()
{
  // Renormalize density function so that it sums to unity. Note that we save
  // the integral of the distribution so that if it is used as part of another
  // distribution (e.g., Mixture), we know its relative strength.
  integral_ = std::accumulate(prob_.begin(), prob_.end(), 0.0);
  for (auto& p_i : prob_) {
    p_i /= integral_;
  }
}

//==============================================================================
// Discrete implementation
//==============================================================================

Discrete::Discrete(pugi::xml_node node) : di_(node)
{
  auto params = get_node_array<double>(node, "parameters");

  std::size_t n = params.size() / 2;

  x_.assign(params.begin(), params.begin() + n);
}

Discrete::Discrete(const double* x, const double* p, size_t n) : di_({p, n})
{

  x_.assign(x, x + n);
}

double Discrete::sample(uint64_t* seed) const
{
  return x_[di_.sample(seed)];
}

//==============================================================================
// Uniform implementation
//==============================================================================

Uniform::Uniform(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 2) {
    fatal_error("Uniform distribution must have two "
                "parameters specified.");
  }

  a_ = params.at(0);
  b_ = params.at(1);
}

double Uniform::sample(uint64_t* seed) const
{
  return a_ + prn(seed) * (b_ - a_);
}

//==============================================================================
// PowerLaw implementation
//==============================================================================

PowerLaw::PowerLaw(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 3) {
    fatal_error("PowerLaw distribution must have three "
                "parameters specified.");
  }

  const double a = params.at(0);
  const double b = params.at(1);
  const double n = params.at(2);

  offset_ = std::pow(a, n + 1);
  span_ = std::pow(b, n + 1) - offset_;
  ninv_ = 1 / (n + 1);
}

double PowerLaw::sample(uint64_t* seed) const
{
  return std::pow(offset_ + prn(seed) * span_, ninv_);
}

//==============================================================================
// Maxwell implementation
//==============================================================================

Maxwell::Maxwell(pugi::xml_node node)
{
  theta_ = std::stod(get_node_value(node, "parameters"));
}

double Maxwell::sample(uint64_t* seed) const
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

double Watt::sample(uint64_t* seed) const
{
  return watt_spectrum(a_, b_, seed);
}

//==============================================================================
// Normal implementation
//==============================================================================
Normal::Normal(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() != 2) {
    openmc::fatal_error("Normal energy distribution must have two "
                        "parameters specified.");
  }

  mean_value_ = params.at(0);
  std_dev_ = params.at(1);
}

double Normal::sample(uint64_t* seed) const
{
  return normal_variate(mean_value_, std_dev_, seed);
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
      openmc::fatal_error(
        "Unsupported interpolation type for distribution: " + temp);
    }
  } else {
    interp_ = Interpolation::histogram;
  }

  // Read and initialize tabular distribution
  auto params = get_node_array<double>(node, "parameters");
  std::size_t n = params.size() / 2;
  const double* x = params.data();
  const double* p = x + n;
  init(x, p, n);
}

Tabular::Tabular(const double* x, const double* p, int n, Interpolation interp,
  const double* c)
  : interp_ {interp}
{
  init(x, p, n, c);
}

void Tabular::init(
  const double* x, const double* p, std::size_t n, const double* c)
{
  // Copy x/p arrays into vectors
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
    std::copy(c, c + n, std::back_inserter(c_));
  } else {
    c_.resize(n);
    c_[0] = 0.0;
    for (int i = 1; i < n; ++i) {
      if (interp_ == Interpolation::histogram) {
        c_[i] = c_[i - 1] + p_[i - 1] * (x_[i] - x_[i - 1]);
      } else if (interp_ == Interpolation::lin_lin) {
        c_[i] = c_[i - 1] + 0.5 * (p_[i - 1] + p_[i]) * (x_[i] - x_[i - 1]);
      }
    }
  }

  // Normalize density and distribution functions. Note that we save the
  // integral of the distribution so that if it is used as part of another
  // distribution (e.g., Mixture), we know its relative strength.
  integral_ = c_[n - 1];
  for (int i = 0; i < n; ++i) {
    p_[i] = p_[i] / integral_;
    c_[i] = c_[i] / integral_;
  }
}

double Tabular::sample(uint64_t* seed) const
{
  // Sample value of CDF
  double c = prn(seed);

  // Find first CDF bin which is above the sampled value
  double c_i = c_[0];
  int i;
  std::size_t n = c_.size();
  for (i = 0; i < n - 1; ++i) {
    if (c <= c_[i + 1])
      break;
    c_i = c_[i + 1];
  }

  // Determine bounding PDF values
  double x_i = x_[i];
  double p_i = p_[i];

  if (interp_ == Interpolation::histogram) {
    // Histogram interpolation
    if (p_i > 0.0) {
      return x_i + (c - c_i) / p_i;
    } else {
      return x_i;
    }
  } else {
    // Linear-linear interpolation
    double x_i1 = x_[i + 1];
    double p_i1 = p_[i + 1];

    double m = (p_i1 - p_i) / (x_i1 - x_i);
    if (m == 0.0) {
      return x_i + (c - c_i) / p_i;
    } else {
      return x_i +
             (std::sqrt(std::max(0.0, p_i * p_i + 2 * m * (c - c_i))) - p_i) /
               m;
    }
  }
}

//==============================================================================
// Equiprobable implementation
//==============================================================================

double Equiprobable::sample(uint64_t* seed) const
{
  std::size_t n = x_.size();

  double r = prn(seed);
  int i = std::floor((n - 1) * r);

  double xl = x_[i];
  double xr = x_[i + i];
  return xl + ((n - 1) * r - i) * (xr - xl);
}

//==============================================================================
// Mixture implementation
//==============================================================================

Mixture::Mixture(pugi::xml_node node)
{
  double cumsum = 0.0;
  for (pugi::xml_node pair : node.children("pair")) {
    // Check that required data exists
    if (!pair.attribute("probability"))
      fatal_error("Mixture pair element does not have probability.");
    if (!pair.child("dist"))
      fatal_error("Mixture pair element does not have a distribution.");

    // cummulative sum of probabilities
    double p = std::stod(pair.attribute("probability").value());

    // Save cummulative probability and distribution
    auto dist = distribution_from_xml(pair.child("dist"));
    cumsum += p * dist->integral();

    distribution_.push_back(std::make_pair(cumsum, std::move(dist)));
  }

  // Save integral of distribution
  integral_ = cumsum;

  // Normalize cummulative probabilities to 1
  for (auto& pair : distribution_) {
    pair.first /= cumsum;
  }
}

double Mixture::sample(uint64_t* seed) const
{
  // Sample value of CDF
  const double p = prn(seed);

  // find matching distribution
  const auto it = std::lower_bound(distribution_.cbegin(), distribution_.cend(),
    p, [](const DistPair& pair, double p) { return pair.first < p; });

  // This should not happen. Catch it
  Ensures(it != distribution_.cend());

  // Sample the chosen distribution
  return it->second->sample(seed);
}

//==============================================================================
// Helper function
//==============================================================================

UPtrDist distribution_from_xml(pugi::xml_node node)
{
  if (!check_for_node(node, "type"))
    openmc::fatal_error("Distribution type must be specified.");

  // Determine type of distribution
  std::string type = get_node_value(node, "type", true, true);

  // Allocate extension of Distribution
  UPtrDist dist;
  if (type == "uniform") {
    dist = UPtrDist {new Uniform(node)};
  } else if (type == "powerlaw") {
    dist = UPtrDist {new PowerLaw(node)};
  } else if (type == "maxwell") {
    dist = UPtrDist {new Maxwell(node)};
  } else if (type == "watt") {
    dist = UPtrDist {new Watt(node)};
  } else if (type == "normal") {
    dist = UPtrDist {new Normal(node)};
  } else if (type == "discrete") {
    dist = UPtrDist {new Discrete(node)};
  } else if (type == "tabular") {
    dist = UPtrDist {new Tabular(node)};
  } else if (type == "mixture") {
    dist = UPtrDist {new Mixture(node)};
  } else if (type == "muir") {
    openmc::fatal_error(
      "'muir' distributions are now specified using the openmc.stats.muir() "
      "function in Python. Please regenerate your XML files.");
  } else {
    openmc::fatal_error("Invalid distribution type: " + type);
  }
  return dist;
}

} // namespace openmc
