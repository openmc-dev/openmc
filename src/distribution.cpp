#include "openmc/distribution.h"

#include <algorithm> // for copy
#include <array>
#include <cmath>     // for sqrt, floor, max
#include <iterator>  // for back_inserter
#include <numeric>   // for accumulate
#include <stdexcept> // for runtime_error
#include <string>    // for string, stod

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

// PDF evaluation not supported for all distribution types
double Distribution::evaluate(double x) const
{
  throw std::runtime_error(
    "PDF evaluation not implemented for this distribution type.");
}

//==============================================================================
// DiscreteIndex implementation
//==============================================================================

DiscreteIndex::DiscreteIndex(pugi::xml_node node)
{
  auto params = get_node_array<double>(node, "parameters");
  std::size_t n = params.size() / 2;

  assign({params.data() + n, n});
}

DiscreteIndex::DiscreteIndex(span<const double> p)
{
  assign(p);
}

void DiscreteIndex::assign(span<const double> p)
{
  prob_.assign(p.begin(), p.end());

  this->init_alias();
  this->init_wgt();
}

void DiscreteIndex::init_alias()
{
  normalize();

  // record user input normalized distribution prob_actual for RR/source bias
  prob_actual_ = prob_;

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

void DiscreteIndex::init_wgt()
{
  wgt_.assign(prob_.size(), 1.0);
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

void DiscreteIndex::apply_bias(span<const double> b)
{
  // Replace the probability vector with that from the bias distribution.
  prob_.assign(b.begin(), b.end());
  if (prob_.size() != prob_actual_.size()) {
    openmc::fatal_error(
      "Size mismatch: Attempted to bias Discrete distribution with " +
      std::to_string(prob_actual_.size()) +
      " probability entries using a "
      "Discrete distribution with " +
      std::to_string(prob_.size()) +
      " entries. Please ensure distributions have the same size.");
  }

  // Normalize biased probability vector and populate weight table.
  normalize();
  for (std::size_t i = 0; i < wgt_.size(); ++i) {
    if (prob_[i] == 0.0) {
      // Allow nonzero entries in original distribution to be given zero
      // sampling probability in the biased distribution.
      wgt_[i] = INFTY;
    } else {
      wgt_[i] = prob_actual_[i] / prob_[i];
    }
  }

  // Reconstruct alias table for sampling from the biased distribution.
  // Values from unbiased prob_actual_ may be recovered using weight table.
  this->init_alias();
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

std::pair<double, double> Discrete::sample(uint64_t* seed) const
{
  size_t sample_index = di_.sample(seed);
  return {x_[sample_index], di_.weight()[sample_index]};
}

double Discrete::evaluate(double x) const
{
  // This function is not called when sampling from a Discrete distribution,
  // even if it is biased. This is because Discrete distributions may only
  // be biased by another Discrete distribution. It is only called when a
  // Discrete distribution is used to bias another kind of distribution.
  for (size_t i = 0; i < x_.size(); ++i) {
    if (std::fabs(x_[i] - x) <= FP_PRECISION) {
      return di_.prob_actual()[i];
    }
  }
  return 0.0;
}

void Discrete::set_bias_discrete(pugi::xml_node node)
{
  // Takes the probability vector from a bias distribution and applies it to
  // the existing DiscreteIndex.
  auto bias_params = get_node_array<double>(node, "bias");

  di_.apply_bias(bias_params);
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

std::pair<double, double> Uniform::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    return {a_ + prn(seed) * (b_ - a_), 1.0};
  }
}

double Uniform::evaluate(double x) const
{
  if (x <= a()) {
    return 0.0;
  } else if (x >= b()) {
    return 0.0;
  } else {
    return 1 / (b() - a());
  }
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

double PowerLaw::evaluate(double x) const
{
  if (x <= a()) {
    return 0.0;
  } else if (x >= b()) {
    return 0.0;
  } else {
    int pwr = n() + 1;
    double norm = pwr / span_;
    return norm * std::pow(std::fabs(x), n());
  }
}

std::pair<double, double> PowerLaw::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    return {std::pow(offset_ + prn(seed) * span_, ninv_), 1.0};
  }
}

//==============================================================================
// Maxwell implementation
//==============================================================================

Maxwell::Maxwell(pugi::xml_node node)
{
  theta_ = std::stod(get_node_value(node, "parameters"));
}

std::pair<double, double> Maxwell::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    return {maxwell_spectrum(theta_, seed), 1.0};
  }
}

double Maxwell::evaluate(double x) const
{
  double c = (2.0 / SQRT_PI) * std::pow(theta_, -1.5);
  return c * std::sqrt(x) * std::exp(-x / theta_);
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

std::pair<double, double> Watt::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    return {watt_spectrum(a_, b_, seed), 1.0};
  }
}

double Watt::evaluate(double x) const
{
  double c = std::exp(-a_ * (b_ / 4.0)) / (std::pow(a_, 2.0) * std::sqrt(b_));
  return c * std::exp(-x / a_) * std::sinh(std::sqrt(b_ * x));
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

std::pair<double, double> Normal::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    return {normal_variate(mean_value_, std_dev_, seed), 1.0};
  }
}

double Normal::evaluate(double x) const
{
  return (1.0 / (std::sqrt(2.0 / PI) * std_dev_)) *
         std::exp(-(std::pow((x - mean_value_), 2.0)) /
                  (2.0 * std::pow(std_dev_, 2.0)));
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

  // Read and initialize tabular distribution. If number of parameters is odd,
  // add an extra zero for the 'p' array.
  auto params = get_node_array<double>(node, "parameters");
  if (params.size() % 2 != 0) {
    params.push_back(0.0);
  }
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

std::pair<double, double> Tabular::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
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
        return {(x_i + (c - c_i) / p_i), 1.0};
      } else {
        return {x_i, 1.0};
      }
    } else {
      // Linear-linear interpolation
      double x_i1 = x_[i + 1];
      double p_i1 = p_[i + 1];

      double m = (p_i1 - p_i) / (x_i1 - x_i);
      if (m == 0.0) {
        return {(x_i + (c - c_i) / p_i), 1.0};
      } else {
        return {
          (x_i +
            (std::sqrt(std::max(0.0, p_i * p_i + 2 * m * (c - c_i))) - p_i) /
              m),
          1.0};
      }
    }
  }
}

double Tabular::evaluate(double x) const
{
  int i;

  if (interp_ == Interpolation::histogram) {
    i = std::upper_bound(x_.begin(), x_.end(), x) - x_.begin() - 1;
    if (i < 0 || i >= static_cast<int>(p_.size())) {
      return 0.0;
    } else {
      return p_[i];
    }
  } else {
    i = std::lower_bound(x_.begin(), x_.end(), x) - x_.begin() - 1;

    if (i < 0 || i >= static_cast<int>(p_.size()) - 1) {
      return 0.0;
    } else {
      double x0 = x_[i];
      double x1 = x_[i + 1];
      double p0 = p_[i];
      double p1 = p_[i + 1];

      double t = (x - x0) / (x1 - x0);
      return (1 - t) * p0 + t * p1;
    }
  }
}

//==============================================================================
// Equiprobable implementation
//==============================================================================

std::pair<double, double> Equiprobable::sample(uint64_t* seed) const
{
  if (bias()) {
    auto [val, wgt] = bias()->sample(seed);
    return {val, this->evaluate(val) / bias()->evaluate(val)};
  } else {
    std::size_t n = x_.size();

    double r = prn(seed);
    int i = std::floor((n - 1) * r);

    double xl = x_[i];
    double xr = x_[i + i];
    return {(xl + ((n - 1) * r - i) * (xr - xl)), 1.0};
  }
}

double Equiprobable::evaluate(double x) const
{
  double x_min = *std::min_element(x_.begin(), x_.end());
  double x_max = *std::max_element(x_.begin(), x_.end());

  if (x < x_min || x > x_max) {
    return 0.0;
  } else {
    return 1.0 / (x_max - x_min);
  }
}

//==============================================================================
// Mixture implementation
//==============================================================================

Mixture::Mixture(pugi::xml_node node)
{
  vector<double> probabilities;

  // First pass: collect distributions and their probabilities
  for (pugi::xml_node pair : node.children("pair")) {
    // Check that required data exists
    if (!pair.attribute("probability"))
      fatal_error("Mixture pair element does not have probability.");
    if (!pair.child("dist"))
      fatal_error("Mixture pair element does not have a distribution.");

    // Get probability and distribution
    double p = std::stod(pair.attribute("probability").value());
    auto dist = distribution_from_xml(pair.child("dist"));

    // Weight probability by the distribution's integral
    double weighted_prob = p * dist->integral();
    probabilities.push_back(weighted_prob);
    distributions_.push_back(std::move(dist));
  }

  // Save sum of weighted probabilities
  integral_ = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);

  // Initialize DiscreteIndex with probability vector, which will normalize
  di_.assign(probabilities);
}

void Mixture::set_bias_mixture(pugi::xml_node node)
{
  // Takes the probability vector from a bias distribution and applies it to
  // the existing DiscreteIndex.
  auto bias_params = get_node_array<double>(node, "bias");

  di_.apply_bias(bias_params);
}

std::pair<double, double> Mixture::sample(uint64_t* seed) const
{
  size_t sample_index = di_.sample(seed);

  // Sample the chosen distribution
  std::pair<double, double> sample_pair =
    distribution_[sample_index]->sample(seed);

  return {sample_pair.first, di_.weight()[sample_index] * sample_pair.second};
}
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

  // Check for biasing distribution
  if (check_for_node(node, "bias")) {
    pugi::xml_node bias_node = node.child("bias");

    if (check_for_node(bias_node, "bias")) {
      openmc::fatal_error(
        "Distribution of type " + type +
        " has a bias distribution with its "
        "own bias distribution. Please ensure bias distributions do not have "
        "their own bias.");
    } else if (type == "discrete") {
      static_cast<Discrete*>(dist.get())->set_bias_discrete(node);
    } else if (type == "mixture") {
      static_cast<Mixture*>(dist.get())->set_bias_mixture(node);
    } else {
      UPtrDist bias = distribution_from_xml(bias_node);
      dist->set_bias(std::move(bias));
    }
  }

  return dist;
}

} // namespace openmc
