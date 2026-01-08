//! \file distribution.h
//! Univariate probability distributions

#ifndef OPENMC_DISTRIBUTION_H
#define OPENMC_DISTRIBUTION_H

#include <cstddef> // for size_t

#include "pugixml.hpp"

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/span.h"
#include "openmc/vector.h" // for vector

namespace openmc {

//==============================================================================
// Helper function for computing importance weights from biased sampling
//==============================================================================

//! Compute importance weights for biased sampling
//! \param p Unnormalized original probability vector
//! \param b Unnormalized bias probability vector
//! \return Vector of importance weights (p_norm[i] / b_norm[i])
vector<double> compute_importance_weights(
  const vector<double>& p, const vector<double>& b);

//==============================================================================
//! Abstract class representing a univariate probability distribution
//==============================================================================

class Distribution {
public:
  virtual ~Distribution() = default;

  //! Sample a value from the distribution, handling biasing automatically
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled value, importance weight)
  virtual std::pair<double, double> sample(uint64_t* seed) const;

  //! Evaluate pdf at a point
  //! \return Value of pdf at a point
  virtual double evaluate(double x) const;

  //! Return integral of distribution
  //! \return Integral of distribution
  virtual double integral() const { return 1.0; };

  //! Set bias distribution
  virtual void set_bias(std::unique_ptr<Distribution> bias)
  {
    bias_ = std::move(bias);
  }

  const Distribution* bias() const { return bias_.get(); }

protected:
  //! Sample without bias handling - each derived class implements core sampling
  //! \param seed Pseudorandom number seed pointer
  //! \return sampled value
  virtual double sample_unbiased(uint64_t* seed) const = 0;

  //! Read bias distribution from XML and call set_bias
  //! \param node XML node that may contain a bias child element
  void read_bias_from_xml(pugi::xml_node node);

  // Biasing distribution
  unique_ptr<Distribution> bias_;
};

using UPtrDist = unique_ptr<Distribution>;

//! Return univariate probability distribution specified in XML file
//! \param[in] node XML node representing distribution
//! \return Unique pointer to distribution
UPtrDist distribution_from_xml(pugi::xml_node node);

//==============================================================================
//! A discrete distribution index (probability mass function)
//==============================================================================

class DiscreteIndex {
public:
  DiscreteIndex() {};
  DiscreteIndex(pugi::xml_node node);
  DiscreteIndex(span<const double> p);

  void assign(span<const double> p);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled index
  size_t sample(uint64_t* seed) const;

  // Properties
  const vector<double>& prob() const { return prob_; }
  const vector<size_t>& alias() const { return alias_; }
  double integral() const { return integral_; }

private:
  vector<double> prob_;  //!< Acceptance probabilities for alias method
  vector<size_t> alias_; //!< Alias table
  double integral_;      //!< Integral of distribution before normalization

  //! Normalize distribution so that probabilities sum to unity
  void normalize();

  //! Initialize alias table for sampling
  void init_alias();
};

//==============================================================================
//! A discrete distribution (probability mass function)
//==============================================================================

class Discrete : public Distribution {
public:
  explicit Discrete(pugi::xml_node node);
  Discrete(const double* x, const double* p, size_t n);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled value, sample weight)
  std::pair<double, double> sample(uint64_t* seed) const override;

  double integral() const override { return di_.integral(); };

  //! Override set_bias as no-op (bias handled in constructor)
  void set_bias(std::unique_ptr<Distribution> bias) override {}

  // Properties
  const vector<double>& x() const { return x_; }
  const vector<double>& prob() const { return di_.prob(); }
  const vector<size_t>& alias() const { return di_.alias(); }
  const vector<double>& weight() const { return weight_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  vector<double> x_;      //!< Possible outcomes
  vector<double> weight_; //!< Importance weights (empty if unbiased)
  DiscreteIndex di_; //!< Discrete probability distribution of outcome indices
};

//==============================================================================
//! Uniform distribution over the interval [a,b]
//==============================================================================

class Uniform : public Distribution {
public:
  explicit Uniform(pugi::xml_node node);
  Uniform(double a, double b) : a_ {a}, b_ {b} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  double a() const { return a_; }
  double b() const { return b_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  double a_; //!< Lower bound of distribution
  double b_; //!< Upper bound of distribution
};

//==============================================================================
//! PowerLaw distribution over the interval [a,b] with exponent n : p(x)=c x^n
//==============================================================================

class PowerLaw : public Distribution {
public:
  explicit PowerLaw(pugi::xml_node node);
  PowerLaw(double a, double b, double n)
    : offset_ {std::pow(a, n + 1)}, span_ {std::pow(b, n + 1) - offset_},
      ninv_ {1 / (n + 1)} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  double a() const { return std::pow(offset_, ninv_); }
  double b() const { return std::pow(offset_ + span_, ninv_); }
  double n() const { return 1 / ninv_ - 1; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  //! Store processed values in object to allow for faster sampling
  double offset_; //!< a^(n+1)
  double span_;   //!< b^(n+1) - a^(n+1)
  double ninv_;   //!< 1/(n+1)
};

//==============================================================================
//! Maxwellian distribution of form c*sqrt(E)*exp(-E/theta)
//==============================================================================

class Maxwell : public Distribution {
public:
  explicit Maxwell(pugi::xml_node node);
  Maxwell(double theta) : theta_ {theta} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  double theta() const { return theta_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  double theta_; //!< Factor in exponential [eV]
};

//==============================================================================
//! Watt fission spectrum with form c*exp(-E/a)*sinh(sqrt(b*E))
//==============================================================================

class Watt : public Distribution {
public:
  explicit Watt(pugi::xml_node node);
  Watt(double a, double b) : a_ {a}, b_ {b} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  double a() const { return a_; }
  double b() const { return b_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  double a_; //!< Factor in exponential [eV]
  double b_; //!< Factor in square root [1/eV]
};

//==============================================================================
//! Normal distributions with form 1/2*std_dev*sqrt(pi) exp
//! (-(e-E0)/2*std_dev)^2
//==============================================================================

class Normal : public Distribution {
public:
  explicit Normal(pugi::xml_node node);
  Normal(double mean_value, double std_dev)
    : mean_value_ {mean_value}, std_dev_ {std_dev} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  double mean_value() const { return mean_value_; }
  double std_dev() const { return std_dev_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  double mean_value_; //!< middle of distribution [eV]
  double std_dev_;    //!< standard deviation [eV]
};

//==============================================================================
//! Histogram or linear-linear interpolated tabular distribution
//==============================================================================

class Tabular : public Distribution {
public:
  explicit Tabular(pugi::xml_node node);
  Tabular(const double* x, const double* p, int n, Interpolation interp,
    const double* c = nullptr);

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  // properties
  vector<double>& x() { return x_; }
  const vector<double>& x() const { return x_; }
  const vector<double>& p() const { return p_; }
  Interpolation interp() const { return interp_; }
  double integral() const override { return integral_; };

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  vector<double> x_;     //!< tabulated independent variable
  vector<double> p_;     //!< tabulated probability density
  vector<double> c_;     //!< cumulative distribution at tabulated values
  Interpolation interp_; //!< interpolation rule
  double integral_;      //!< Integral of distribution

  //! Initialize tabulated probability density function
  //! \param x Array of values for independent variable
  //! \param p Array of tabulated probabilities
  //! \param n Number of tabulated values
  void init(
    const double* x, const double* p, std::size_t n, const double* c = nullptr);
};

//==============================================================================
//! Equiprobable distribution
//==============================================================================

class Equiprobable : public Distribution {
public:
  explicit Equiprobable(pugi::xml_node node);
  Equiprobable(const double* x, int n) : x_ {x, x + n} {};

  //! Evaluate probability density at a point
  //! \param x Point to evaluate f(x)
  //! \return f(x)
  double evaluate(double x) const override;

  const vector<double>& x() const { return x_; }

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  vector<double> x_; //! Possible outcomes
};

//==============================================================================
//! Mixture distribution
//==============================================================================

class Mixture : public Distribution {
public:
  explicit Mixture(pugi::xml_node node);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return (sampled value, sample weight)
  std::pair<double, double> sample(uint64_t* seed) const override;

  double integral() const override { return integral_; }

  //! Override set_bias as no-op (bias handled in constructor)
  void set_bias(std::unique_ptr<Distribution> bias) override {}

protected:
  double sample_unbiased(uint64_t* seed) const override;

private:
  vector<UPtrDist> distribution_; //!< Sub-distributions
  vector<double> weight_; //!< Importance weights for component selection
  DiscreteIndex di_;      //!< Discrete probability distribution of indices
  double integral_;       //!< Integral of distribution
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_H
