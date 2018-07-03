#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>

#include "constants.h"

namespace openmc {

class Distribution {
public:
  virtual double sample() = 0;
  virtual ~Distribution() = default;
};


class Discrete : public Distribution {
public:
  Discrete(const double* x, const double* p, int n);
  double sample();
private:
  std::vector<double> x_;
  std::vector<double> p_;
};


class Uniform : public Distribution {
public:
  Uniform(double a, double b) : a_{a}, b_{b} {};
  double sample();
private:
  double a_;
  double b_;
};


class Maxwell : public Distribution {
public:
  Maxwell(double theta) : theta_{theta} { };
  double sample();
private:
  double theta_;
};


class Watt : public Distribution {
public:
  Watt(double a, double b) : a_{a}, b_{b} { };
  double sample();
private:
  double a_;
  double b_;
};


class Tabular : public Distribution {
public:
  Tabular(const double* x, const double* p, int n, Interpolation interp);
  double sample();
private:
  std::vector<double> x_;
  std::vector<double> p_;
  std::vector<double> c_;
  Interpolation interp_;
};


class Equiprobable : public Distribution {
public:
  Equiprobable(const double* x, int n) : x_{x, x+n} { };
  double sample();
private:
  std::vector<double> x_;
};

} // namespace openmc

#endif // DISTRIBUTION_H
