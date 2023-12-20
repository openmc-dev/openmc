//! \file distribution.h
//! Univariate probability distributions

#ifndef OPENMC_DISTRIBUTION_H
#define OPENMC_DISTRIBUTION_H

#include <cstddef> // for size_t

#include "pugixml.hpp"
#include <gsl/gsl-lite.hpp>

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/vector.h" // for vector

namespace openmc {

//==============================================================================
//! Abstract class representing a univariate probability distribution
//==============================================================================

class Distribution {
public:
  virtual ~Distribution() = default;
  virtual double sample(uint64_t* seed) const = 0;

  //! Return integral of distribution
  //! \return Integral of distribution
  virtual double integral() const { return 1.0; };
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
  DiscreteIndex(gsl::span<const double> p);

  void assign(gsl::span<const double> p);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  size_t sample(uint64_t* seed) const;

  // Properties
  const vector<double>& prob() const { return prob_; }
  const vector<size_t>& alias() const { return alias_; }
  double integral() const { return integral_; }

private:
  vector<double> prob_; //!< Probability of accepting the uniformly sampled bin,
                        //!< mapped to alias method table
  vector<size_t> alias_; //!< Alias table
  double integral_;      //!< Integral of distribution

  //! Normalize distribution so that probabilities sum to unity
  void normalize();

  //! Initialize alias tables for distribution
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
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double integral() const override { return di_.integral(); };

  // Properties
  const vector<double>& x() const { return x_; }
  const vector<double>& prob() const { return di_.prob(); }
  const vector<size_t>& alias() const { return di_.alias(); }

private:
  vector<double> x_; //!< Possible outcomes
  DiscreteIndex di_; //!< discrete probability distribution of
                     //!< outcome indices
};

//==============================================================================
//! Uniform distribution over the interval [a,b]
//==============================================================================

class Uniform : public Distribution {
public:
  explicit Uniform(pugi::xml_node node);
  Uniform(double a, double b) : a_ {a}, b_ {b} {};

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double a() const { return a_; }
  double b() const { return b_; }

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double a() const { return std::pow(offset_, ninv_); }
  double b() const { return std::pow(offset_ + span_, ninv_); }
  double n() const { return 1 / ninv_ - 1; }

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double theta() const { return theta_; }

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double a() const { return a_; }
  double b() const { return b_; }

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double mean_value() const { return mean_value_; }
  double std_dev() const { return std_dev_; }

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  // properties
  vector<double>& x() { return x_; }
  const vector<double>& x() const { return x_; }
  const vector<double>& p() const { return p_; }
  Interpolation interp() const { return interp_; }
  double integral() const override { return integral_; };

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

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  const vector<double>& x() const { return x_; }

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
  //! \return Sampled value
  double sample(uint64_t* seed) const override;

  double integral() const override { return integral_; }

private:
  // Storrage for probability + distribution
  using DistPair = std::pair<double, UPtrDist>;

  vector<DistPair>
    distribution_;  //!< sub-distributions + cummulative probabilities
  double integral_; //!< integral of distribution
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_H
