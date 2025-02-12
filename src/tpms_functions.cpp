#include "openmc/tpms_functions.h"

#include <iostream>

#include "openmc/error.h"

namespace openmc {

// *****************************************************************************
// *   FONCTION FOR TPMS GENERAL DEFINITION
// *****************************************************************************

FunctionForTPMS::FunctionForTPMS(double _xmin, double _xmax, double _ymin,
  double _ymax, double _zmin, double _zmax)
  : xmin(_xmin), xmax(_xmax), ymin(_ymin), ymax(_ymax), zmin(_zmin), zmax(_zmax)
{
  minimalValue = this->get_minimum();
}

// *****************************************************************************
// *   HELPER FUNCTIONS TO INTERPOLATE
// *****************************************************************************

void interpolate_1d(double coord, const std::vector<double>& grid_coords,
  int& lower_index, int& upper_index, double& lower_weight,
  double& upper_weight)
{
  int num_points = grid_coords.size();

  auto it = std::upper_bound(grid_coords.begin(), grid_coords.end(), coord);
  upper_index = std::min(
    static_cast<int>(std::distance(grid_coords.begin(), it)), num_points - 1);
  lower_index = std::max(upper_index - 1, 0);

  double lower_coord = grid_coords[lower_index];
  double upper_coord = grid_coords[upper_index];

  if (coord <= grid_coords.front()) {
    lower_index = upper_index = 0;
    lower_weight = 1.0;
    upper_weight = 0.0;
  } else if (coord >= grid_coords.back()) {
    lower_index = upper_index = num_points - 1;
    lower_weight = 0.0;
    upper_weight = 1.0;
  } else {
    double range = upper_coord - lower_coord;
    lower_weight = (upper_coord - coord) / range;
    upper_weight = (coord - lower_coord) / range;
  }
}

// *****************************************************************************
// *   INTERPOLATION FOR TPMS DEFINITION
// *****************************************************************************

InterpolationForTPMS::InterpolationForTPMS(std::vector<double> _x_grid,
  std::vector<double> _y_grid, std::vector<double> _z_grid,
  xt::xarray<double> _matrix)
  : FunctionForTPMS(_x_grid[0], _x_grid[_x_grid.size() - 1], _y_grid[0],
      _y_grid[_y_grid.size() - 1], _z_grid[0], _z_grid[_z_grid.size() - 1]),
    x_grid(_x_grid), y_grid(_y_grid), z_grid(_z_grid), matrix(_matrix)
{
  minimalValue = xt::amin(matrix)();
  useFirstDerivatives = false;
  useSecondDerivatives = false;
}

double InterpolationForTPMS::fxyz(double x, double y, double z) const
{
  int ix0, ix1, iy0, iy1, iz0, iz1;
  double wx0, wx1, wy0, wy1, wz0, wz1;

  interpolate_1d(x, x_grid, ix0, ix1, wx0, wx1);
  interpolate_1d(y, y_grid, iy0, iy1, wy0, wy1);
  interpolate_1d(z, z_grid, iz0, iz1, wz0, wz1);

  double c00 = matrix(iz0,iy0,ix0) * wx0 + matrix(iz0,iy0,ix1) * wx1;
  double c01 = matrix(iz1,iy0,ix0) * wx0 + matrix(iz1,iy0,ix1) * wx1;
  double c10 = matrix(iz0,iy1,ix0) * wx0 + matrix(iz0,iy1,ix1) * wx1;
  double c11 = matrix(iz1,iy1,ix0) * wx0 + matrix(iz1,iy1,ix1) * wx1;

  double c0 = c00 * wy0 + c10 * wy1;
  double c1 = c01 * wy0 + c11 * wy1;
  return c0 * wz0 + c1 * wz1;
}

} // namespace openmc