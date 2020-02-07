//! \file distribution.h
//! Univariate probability distributions

#ifndef OPENMC_DISTRIBUTION_H
#define OPENMC_DISTRIBUTION_H

#include <cstddef> // for size_t
#include <memory> // for unique_ptr
#include <vector> // for vector

#include "pugixml.hpp"

#include "openmc/constants.h"

namespace openmc {

//==============================================================================
//! Abstract class representing a univariate probability distribution
//==============================================================================

class Distribution {
public:
  virtual ~Distribution() = default;
  virtual double sample(uint64_t* seed) const = 0;
};

//==============================================================================
//! A discrete distribution (probability mass function)
//==============================================================================

class Discrete : public Distribution {
public:
  explicit Discrete(pugi::xml_node node);
  Discrete(const double* x, const double* p, int n);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  // Properties
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& p() const { return p_; }
private:
  std::vector<double> x_; //!< Possible outcomes
  std::vector<double> p_; //!< Probability of each outcome

  //! Normalize distribution so that probabilities sum to unity
  void normalize();
};

//==============================================================================
//! Uniform distribution over the interval [a,b]
//==============================================================================

class Uniform : public Distribution {
public:
  explicit Uniform(pugi::xml_node node);
  Uniform(double a, double b) : a_{a}, b_{b} {};

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  double a() const { return a_; }
  double b() const { return b_; }
private:
  double a_; //!< Lower bound of distribution
  double b_; //!< Upper bound of distribution
};

//==============================================================================
//! Maxwellian distribution of form c*E*exp(-E/theta)
//==============================================================================

class Maxwell : public Distribution {
public:
  explicit Maxwell(pugi::xml_node node);
  Maxwell(double theta) : theta_{theta} { };

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

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
  Watt(double a, double b) : a_{a}, b_{b} { };

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  double a() const { return a_; }
  double b() const { return b_; }
private:
  double a_; //!< Factor in exponential [eV]
  double b_; //!< Factor in square root [1/eV]
};

//==============================================================================
//! Normal distributions with form 1/2*std_dev*sqrt(pi) exp (-(e-E0)/2*std_dev)^2
//==============================================================================

class Normal : public Distribution {
public:
  explicit Normal(pugi::xml_node node);
  Normal(double mean_value, double std_dev) : mean_value_{mean_value}, std_dev_{std_dev} { };

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  double mean_value() const { return mean_value_; }
  double std_dev() const { return std_dev_; }
private:
  double mean_value_;    //!< middle of distribution [eV]
  double std_dev_; //!< standard deviation [eV]
};

//==============================================================================
//! Muir (fusion) spectrum derived from Normal with extra params e0 is mean
//! std dev is sqrt(4*e0*kt/m)
//==============================================================================

class Muir : public Distribution {
public:
  explicit Muir(pugi::xml_node node);
  Muir(double e0, double m_rat, double kt) : e0_{e0}, m_rat_{m_rat}, kt_{kt} { };

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  double e0() const { return e0_; }
  double m_rat() const { return m_rat_; }
  double kt() const { return kt_; }
private:
  // example DT fusion m_rat = 5 (D = 2 + T = 3)
  // ion temp = 20000 eV
  // mean neutron energy 14.08e6 eV
  double e0_;    //!< mean neutron energy [eV]
  double m_rat_; //!< ratio of reactant masses relative to atomic mass unit
  double kt_;    //!< ion temperature [eV]
};

//==============================================================================
//! Histogram or linear-linear interpolated tabular distribution
//==============================================================================

class Tabular : public Distribution {
public:
  explicit Tabular(pugi::xml_node node);
  Tabular(const double* x, const double* p, int n, Interpolation interp,
          const double* c=nullptr);

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  // x property
  std::vector<double>& x() { return x_; }
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& p() const { return p_; }
  Interpolation interp() const { return interp_; }
private:
  std::vector<double> x_; //!< tabulated independent variable
  std::vector<double> p_; //!< tabulated probability density
  std::vector<double> c_; //!< cumulative distribution at tabulated values
  Interpolation interp_;  //!< interpolation rule

  //! Initialize tabulated probability density function
  //! \param x Array of values for independent variable
  //! \param p Array of tabulated probabilities
  //! \param n Number of tabulated values
  void init(const double* x, const double* p, std::size_t n,
            const double* c=nullptr);
};

//==============================================================================
//! Equiprobable distribution
//==============================================================================

class Equiprobable : public Distribution {
public:
  explicit Equiprobable(pugi::xml_node node);
  Equiprobable(const double* x, int n) : x_{x, x+n} { };

  //! Sample a value from the distribution
  //! \param seed Pseudorandom number seed pointer
  //! \return Sampled value
  double sample(uint64_t* seed) const;

  const std::vector<double>& x() const { return x_; }
private:
  std::vector<double> x_; //! Possible outcomes
};


using UPtrDist = std::unique_ptr<Distribution>;

//! Return univariate probability distribution specified in XML file
//! \param[in] node XML node representing distribution
//! \return Unique pointer to distribution
UPtrDist distribution_from_xml(pugi::xml_node node);

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_H
