#include "openmc/implicit_solvers.h"

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"

namespace openmc {

double NaiveLipschitz::solve(const Implicit& function, Position r, Direction u,
  double t0, double t1, double isovalue) const
{
  const double L = function.compute_lipschitz(r, u, t0, t1);

  // Constant function — root only if already within tolerance at t0.
  if (L == 0.0) {
    double f0 = function.along_ray(t0, r, u) - isovalue;
    return (std::abs(f0) < atol_) ? t0 : INFTY;
  }

  double t_curr = t0;
  double f_curr = function.along_ray(t_curr, r, u) - isovalue;

  for (int i = 0; i < max_iter_; ++i) {
    // Safe Lipschitz step: f can only change by L*step over this distance,
    // so no root is possible before t_curr + |f_curr| / L.
    const double step = std::abs(f_curr) / L;
    const double t_next = t_curr + step;

    // Stepped past the end of the interval — check endpoint.
    if (t_next >= t1) {
      double f1 = function.along_ray(t1, r, u) - isovalue;
      if (f_curr * f1 < 0.0)
        return bisect(function, r, u, t_curr, t1, f_curr, f1, isovalue);
      return INFTY;
    }

    const double f_next = function.along_ray(t_next, r, u) - isovalue;

    if (std::abs(f_next) < L * atol_)
      return t_next;

    // Sign change detected — root is bracketed in [t_curr, t_next].
    if (f_curr * f_next < 0.0)
      return bisect(function, r, u, t_curr, t_next, f_curr, f_next, isovalue);

    t_curr = t_next;
    f_curr = f_next;
  }

  // Max iterations reached — skip this surface.
  warning(fmt::format("NaiveLipschitz reached max iterations ({})."
                      " The surface will be skipped. Results could be biased.",
    max_iter_));
  return INFTY;
}

double NaiveLipschitz::bisect(const Implicit& function, Position r, Direction u,
  double ta, double tb, double fa, double fb, double isovalue) const
{
  // Invariant: fa and fb have opposite signs, root is in (ta, tb).
  // tb is always on the destination side of the crossing.
  // We return tb rather than the midpoint so that the returned position
  // is unambiguously on the destination side — this prevents sense()
  // from placing the particle back in the cell it just left.
  while ((tb - ta) > atol_) {
    const double tm = 0.5 * (ta + tb);
    const double fm = function.along_ray(tm, r, u) - isovalue;
    if (fa * fm < 0.0) {
      tb = tm;
      fb = fm;
    } else {
      ta = tm;
      fa = fm;
    }
  }
  return tb;
}

} // namespace openmc