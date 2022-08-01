
#ifndef OPENMC_INTERPOLATE_H
#define OPENMC_INTERPOLATE_H

#include <cmath>
#include <vector>

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
  return y0 + log(x / x0) / log(x1 / x0) * (y1 - y0);
}

inline double interpolate_log_lin(
  double x0, double x1, double y0, double y1, double x)
{
  return y0 * exp((x - x0) / (x1 - x0) * log(y1 / y0));
}

inline double interpolate_log_log(
  double x0, double x1, double y0, double y1, double x)
{
  double f = log(x / x0) / log(x1 / x0);
  return y0 * exp(f * log(y1 / y0));
}

inline double interpolate(const std::vector<double>& xs,
  const std::vector<double>& ys, int idx, double x,
  Interpolation i = Interpolation::lin_lin)
{
  switch (i) {
  case Interpolation::lin_lin:
    return interpolate_lin_lin(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::log_log:
    return interpolate_log_log(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::lin_log:
    return interpolate_lin_log(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  case Interpolation::log_lin:
    return interpolate_log_lin(xs[idx], xs[idx + 1], ys[idx], ys[idx + 1], x);
  default:
    fatal_error("Unsupported interpolation");
  }
}

inline double interpolate(const std::vector<double>& xs,
  const std::vector<double>& ys, double x,
  Interpolation i = Interpolation::lin_lin)
{
  int idx = lower_bound_index(xs.begin(), xs.end(), x);
  return interpolate(xs, ys, idx, x, i);
}

} // namespace openmc

#endif