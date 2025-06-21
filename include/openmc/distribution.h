//! \file distribution.h
//! Univariate probability distributions

#ifndef OPENMC_DISTRIBUTION_H
#define OPENMC_DISTRIBUTION_H

#include <algorithm>
#include <cstddef> // for size_t

#include "pugixml.hpp"

#include "openmc/constants.h"
#include "openmc/math_functions.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/random_lcg.h"
#include "openmc/span.h"
#include "openmc/vector.h" // for vector

namespace openmc {

//==============================================================================
//! Abstract class representing a univariate probability distribution
//==============================================================================

class Distribution {
public:
  virtual ~Distribution() = default;

  virtual double sample(uint64_t* seed) const = 0;

  //! Return integral of distribution over whole range
  //! \return Integral of distribution over whole range
  virtual double integral() const { return 1.0; };

  //! Return integral of distribution over finite interval
  //! \return Integral of distribution over finite interval
  double integral(double x0, double x1) const;

  struct Support {
    double first;
    double second;
  };

  //! Return range of possible values
  virtual Support support() const = 0;

  int32_t dims() const { return dims_; }

protected:
  int32_t dims_ {
    C_NONE}; //! Number of uniform random samples needed for sampling
};

class FixedDistribution : public Distribution {
public:
  int32_t dims() const { return dims_; }
  virtual double sample(vector<double>::iterator x) const = 0;
  double sample(uint64_t* seed) const override
  {
    vector<double> x;
    for (auto i = 0; i < dims(); i++) {
      x.emplace_back(prn(seed));
    }
    return this->sample(begin(x));
  }
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
  //! \return Sampled value
  size_t sample(uint64_t* seed) const
  {
    size_t n = prob_.size();
    if (n > 1) {
      vector<double> x = {prn(seed), prn(seed)};
      return sample(begin(x));
    } else {
      return 0;
    }
  }
  size_t sample(vector<double>::iterator x) const;

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

class Discrete : public FixedDistribution {
public:
  explicit Discrete(pugi::xml_node node);
  Discrete(const double* x, const double* p, size_t n);

  using FixedDistribution::sample;
  double sample(vector<double>::iterator x) const override;

  double integral() const override { return di_.integral(); };

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = *std::min_element(begin(x_), end(x_));
    sup.second = *std::max_element(begin(x_), end(x_));
    return sup;
  }

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

class Uniform : public FixedDistribution {
public:
  explicit Uniform(pugi::xml_node node);
  Uniform(double a, double b) : a_ {a}, b_ {b} { dims_ = 1; };

  using FixedDistribution::sample;
  double sample(vector<double>::iterator x) const override;

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = a_;
    sup.second = b_;
    return sup;
  }

  double a() const { return a_; }
  double b() const { return b_; }

private:
  double a_; //!< Lower bound of distribution
  double b_; //!< Upper bound of distribution
};

//==============================================================================
//! PowerLaw distribution over the interval [a,b] with exponent n : p(x)=c x^n
//==============================================================================

class PowerLaw : public FixedDistribution {
public:
  explicit PowerLaw(pugi::xml_node node);
  PowerLaw(double a, double b, double n)
    : offset_ {std::pow(a, n + 1)}, span_ {std::pow(b, n + 1) - offset_},
      ninv_ {1 / (n + 1)}
  {
    dims_ = 1;
  };

  using FixedDistribution::sample;
  double sample(vector<double>::iterator x) const override;

  double integral(double x0, double x1) const;

  double a() const { return std::pow(offset_, ninv_); }
  double b() const { return std::pow(offset_ + span_, ninv_); }
  double n() const { return 1 / ninv_ - 1; }

  Support support() const override
  {
    Support sup;
    sup.first = a();
    sup.second = b();
    return sup;
  }

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

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = 0.0;
    sup.second = INFTY;
    return sup;
  }

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

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = 0.0;
    sup.second = INFTY;
    return sup;
  }

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

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = 0.0;
    sup.second = INFTY;
    return sup;
  }

  double mean_value() const { return mean_value_; }
  double std_dev() const { return std_dev_; }

private:
  double mean_value_; //!< middle of distribution [eV]
  double std_dev_;    //!< standard deviation [eV]
};

//==============================================================================
//! Histogram or linear-linear interpolated tabular distribution
//==============================================================================

class Tabular : public FixedDistribution {
public:
  explicit Tabular(pugi::xml_node node);
  Tabular(const double* x, const double* p, int n, Interpolation interp,
    const double* c = nullptr);

  using FixedDistribution::sample;
  double sample(vector<double>::iterator x) const override;

  // properties
  vector<double>& x() { return x_; }
  const vector<double>& x() const { return x_; }
  const vector<double>& p() const { return p_; }
  Interpolation interp() const { return interp_; }

  double integral() const override { return integral_; }

  Support support() const override
  {
    Support sup;
    sup.first = *std::min_element(begin(x_), end(x_));
    sup.second = *std::max_element(begin(x_), end(x_));
    return sup;
  }

  double integral(double x0, double x1) const;

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

class Equiprobable : public FixedDistribution {
public:
  explicit Equiprobable(pugi::xml_node node);
  Equiprobable(const double* x, int n) : x_ {x, x + n} { dims_ = 1; };

  using FixedDistribution::sample;
  double sample(vector<double>::iterator x) const override;

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = *std::min_element(begin(x_), end(x_));
    sup.second = *std::max_element(begin(x_), end(x_));
    return sup;
  }

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

  double integral(double x0, double x1) const;

  Support support() const override
  {
    Support sup;
    sup.first = INFTY;
    sup.second = -INFTY;
    for (auto& item : distribution_) {
      Support dsup = item.second->support();
      sup.first = std::min(sup.first, dsup.first);
      sup.second = std::max(sup.second, dsup.second);
    }
    return sup;
  }

private:
  // Storage for probability + distribution
  using DistPair = std::pair<double, UPtrDist>;

  vector<DistPair>
    distribution_;  //!< sub-distributions + cummulative probabilities
  double integral_; //!< integral of distribution
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_H
