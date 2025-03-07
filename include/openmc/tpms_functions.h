#ifndef TPMS_FUNCTIONS
#define TPMS_FUNCTIONS

#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"

#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Generic Function class to defined pitch and isovalue as function of x,y
//! and z. 
//==============================================================================

class FunctionForTPMS {
public:
  FunctionForTPMS(double _xmin, double _xmax, double _ymin, double _ymax,
    double _zmin, double _zmax);
  bool is_first_derivatives() { return useFirstDerivatives; };

public:
  virtual double fxyz(double x, double y, double z) const { return 0.; };
  virtual std::vector<double> first_partial_derivatives(
    double x, double y, double z) const
  {
    return {0., 0., 0.};
  };
  virtual std::vector<std::vector<double>> second_partial_derivatives(
    double x, double y, double z) const
  {
    return {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
  };
  virtual double get_minimum() const { return 0.; };

public:
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double minimalValue;
  bool useFirstDerivatives = false;
  bool useSecondDerivatives = false;
};

//==============================================================================
//! Interpolation Function class to interpolate values of pitch and isovalue
//! from a 3D matrix with given of x,y and z grids.
//==============================================================================

class InterpolationForTPMS : public FunctionForTPMS {
public:
  InterpolationForTPMS(std::vector<double> _x_grid, std::vector<double> _y_grid,
    std::vector<double> _z_grid, xt::xarray<double> _matrix);
  double fxyz(double x, double y, double z) const;

public:
  std::vector<double> x_grid;
  std::vector<double> y_grid;
  std::vector<double> z_grid;
  xt::xarray<double> matrix;
};

} // namespace openmc

#endif