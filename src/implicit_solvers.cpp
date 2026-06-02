#include "openmc/implicit_solvers.h"

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/error.h"

namespace openmc {

std::unique_ptr<ImplicitSolver> ImplicitSolver::create(
  const std::string& name, int max_iter, double atol, double ftol)
{
  if (name == "naive")
    return std::make_unique<NaiveLipschitz>(atol, ftol, max_iter);
  if (name == "fast")
    return std::make_unique<FastLipschitz>(atol, ftol, max_iter);

  fatal_error("make_implicit_solver: unknown solver name '" + name +
              "'. "
              "Valid options are: 'naive', 'simple', 'fast'.");
  return nullptr; // unreachable — suppresses compiler warning
}

double NaiveLipschitz::solve(const Implicit& function, Position r, Direction u,
  double t0, double t1, double isovalue) const
{
  const double L = function.compute_lipschitz(r, u, t0, t1);

  // Constant function — root only if already within tolerance at t0.
  if (L == 0.0) {
    double f0 = function.along_ray(t0, r, u) - isovalue;
    return (std::abs(f0) < ftol_) ? t0 : INFTY;
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

    if (std::abs(f_next) < ftol_)
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

double FastLipschitz::solve(const Implicit& function, Position r, Direction u,
  double t0, double t1, double isovalue) const
{
  // Stack of intervals to explore.  Stored as (t0, t1) pairs.
  // Because the stack is LIFO, pushing the right half first ensures
  // the left half (smaller t) is always explored first, so the
  // solver finds the smallest root.
  std::vector<std::pair<double, double>> stack;
  stack.reserve(64);
  stack.push_back({t0, t1});

  int n_iter = 0;

  while (!stack.empty()) {
    auto [t0, t1] = stack.back();
    stack.pop_back();

    // Interval has collapsed to within tolerance — root is at t0.
    if (t1 - t0 < atol_)
      return t1;

    // Compute tight Lipschitz bound and function range over [t0, t1].
    const double L = function.compute_lipschitz(r, u, t0, t1);

    if (L == 0.0) {
      // Constant on this interval — root only if already on surface.
      double f0 = function.along_ray(t0, r, u) - isovalue;
      if (std::abs(f0) < atol_)
        return t0;
      continue;
    }

    auto [fmin_raw, fmax_raw] = function.compute_f_min_max(r, u, t0, t1);
    const double fmin = fmin_raw - isovalue;
    const double fmax = fmax_raw - isovalue;

    // Both bounds have the same sign — no root in this interval.
    if (fmin * fmax > 0.0)
      continue;

    // Endpoint values (needed to compute safe sub-interval).
    const double f0 = function.along_ray(t0, r, u) - isovalue;
    const double f1 = function.along_ray(t1, r, u) - isovalue;

    // Safe start: root cannot exist before tSafeStart.
    // Derivation: |f(t)| >= |f0| - L*(t-t0), so f cannot cross zero
    // until |f0| - L*(t-t0) <= 0  =>  t >= t0 + |f0|/L.
    // Subtract L*atol_ to account for the tolerance band.
    const double tSafeStart = t0 + std::max(0.0, std::abs(f0) - L * atol_) / L;
    if (tSafeStart - t0 < atol_)
      return tSafeStart; // root is within atol of t0

    // Safe end: root cannot exist after tSafeEnd (symmetric argument).
    const double tSafeEnd = t1 - std::max(0.0, std::abs(f1) - L * atol_) / L;
    if (t1 - tSafeEnd < atol_)
      return t1; // root is within atol of t1

    // Safe interval has collapsed — root is at tSafeStart.
    if (tSafeEnd - tSafeStart < atol_)
      return tSafeEnd;

    if (++n_iter >= max_iter_) {
      warning(
        fmt::format("FastLipschitz reached max iterations ({})."
                    " The surface will be skipped. Results could be biased.",
          max_iter_));
      return INFTY;
    }

    // Split the safe interval and push both halves.
    // Right is pushed first so that left is popped and explored first.
    const double tSplit = 0.5 * (tSafeStart + tSafeEnd);
    stack.push_back({tSplit, tSafeEnd});   // right — explored second
    stack.push_back({tSafeStart, tSplit}); // left  — explored first
  }

  return INFTY;
}

} // namespace openmc