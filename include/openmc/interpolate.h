#ifndef OPENMC_INTERPOLATE_H
#define OPENMC_INTERPOLATE_H

#include <cmath>
#include <vector>

#include <gsl/gsl-lite.hpp>

#include "openmc/error.h"
#include "openmc/search.h"

namespace openmc {

inline double interpolate_lin_lin(
  double x0, double x1, double y0, double y1, double x)
{
  return y0 + (x - x0) / (x1 - x0) * (y1 - y0);
}

inline double interpolate_lin_log(
  double x0, double x1, double y0, double y1, double x)
{
  return y0 + std::log(x / x0) / std::log(x1 / x0) * (y1 - y0);
}

inline double interpolate_log_lin(
  double x0, double x1, double y0, double y1, double x)
{
  return y0 * std::exp((x - x0) / (x1 - x0) * std::log(y1 / y0));
}

inline double interpolate_log_log(
  double x0, double x1, double y0, double y1, double x)
{
  double f = std::log(x / x0) / std::log(x1 / x0);
  return y0 * std::exp(f * std::log(y1 / y0));
}

inline double interpolate_lagrangian(gsl::span<const double> xs,
  gsl::span<const double> ys, int idx, double x, int order)
{
  double output {0.0};

  for (int i = 0; i < order + 1; i++) {
    double numerator {1.0};
    double denominator {1.0};
    for (int j = 0; j < order + 1; j++) {
      if (i == j)
        continue;
      numerator *= (x - xs[idx + j]);
      denominator *= (xs[idx + i] - xs[idx + j]);
    }
    output += (numerator / denominator) * ys[idx + i];
  }

  return output;
}

inline double interpolate(gsl::span<const double> xs,
  gsl::span<const double> ys, double x,
  Interpolation i = Interpolation::lin_lin)
{
  int idx = lower_bound_index(xs.begin(), xs.end(), x);

  if (idx == xs.size())
    idx--;

  switch (i) {
  case Interpolation::histogram:
    return ys[idx];
  case Interpolation::lin_lin:
    return interpolate_lin_lin(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::log_log:
    return interpolate_log_log(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::lin_log:
    return interpolate_lin_log(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::log_lin:
    return interpolate_log_lin(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::quadratic:
    // move back one point if x is in the last interval of the x-grid
    if (idx == xs.size() - 2 && idx > 0)
      idx--;
    return interpolate_lagrangian(xs, ys, idx, x, 2);
  case Interpolation::cubic:
    // if x is not in the first interval of the x-grid, move back one
    if (idx > 0)
      idx--;
    // if the index was the last interval of the x-grid, move it back one more
    if (idx == xs.size() - 3)
      idx--;
    return interpolate_lagrangian(xs, ys, idx, x, 3);
  default:
    fatal_error("Unsupported interpolation");
  }
}

} // namespace openmc

#endif
