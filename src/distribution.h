#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <cstddef> // for size_t
#include <memory> // for unique_ptr
#include <vector> // for vector

#include "constants.h"
#include "pugixml.hpp"

namespace openmc {

//==============================================================================
//! Abstract class representing a univariate probability distribution
//==============================================================================

class Distribution {
public:
  virtual double sample() const = 0;
  virtual ~Distribution() = default;
};

//==============================================================================
//! A discrete distribution (probability mass function)
//==============================================================================

class Discrete : public Distribution {
public:
  explicit Discrete(pugi::xml_node node);
  Discrete(const double* x, const double* p, int n);

  //! Sample a value from the distribution
  double sample() const;
private:
  std::vector<double> x_;
  std::vector<double> p_;

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
  double sample() const;
private:
  double a_;
  double b_;
};

//==============================================================================
//! Maxwellian distribution of form c*E*exp(-E/a)
//==============================================================================

class Maxwell : public Distribution {
public:
  explicit Maxwell(pugi::xml_node node);
  Maxwell(double theta) : theta_{theta} { };

  //! Sample a value from the distribution
  double sample() const;
private:
  double theta_;
};

//==============================================================================
//! Watt fission spectrum with form c*exp(-E/a)*sinh(sqrt(b*E))
//==============================================================================

class Watt : public Distribution {
public:
  explicit Watt(pugi::xml_node node);
  Watt(double a, double b) : a_{a}, b_{b} { };

  //! Sample a value from the distribution
  double sample() const;
private:
  double a_;
  double b_;
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
  double sample() const;
private:
  std::vector<double> x_; //!< tabulated independent variable
  std::vector<double> p_; //!< tabulated probability density
  std::vector<double> c_; //!< cumulative distribution at tabulated values
  Interpolation interp_;  //!< interpolation rule

  //! Initialize tabulated probability density function
  //! @param x Array of values for independent variable
  //! @param p Array of tabulated probabilities
  //! @param n Number of tabulated values
  void init(const double* x, const double* p, std::size_t n,
            const double* c=nullptr);
};

//==============================================================================
//!
//==============================================================================

class Equiprobable : public Distribution {
public:
  explicit Equiprobable(pugi::xml_node node);
  Equiprobable(const double* x, int n) : x_{x, x+n} { };

  //! Sample a value from the distribution
  double sample() const;
private:
  std::vector<double> x_;
};


using UPtrDist = std::unique_ptr<Distribution>;

UPtrDist distribution_from_xml(pugi::xml_node node);

} // namespace openmc

#endif // DISTRIBUTION_H
