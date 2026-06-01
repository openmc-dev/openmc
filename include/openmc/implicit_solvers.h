#ifndef OPENMC_IMPLICIT_SOLVER_H
#define OPENMC_IMPLICIT_SOLVER_H

#include "openmc/implicit.h"
#include "openmc/position.h"

namespace openmc {

//==============================================================================
// ImplicitSolver
//
// Abstract base class for ray-implicit surface intersection solvers.
//
// Contract:
//   - func has already been shifted by isovalue (SurfaceImplicit::evaluate
//     returns f(r) - c, so the solver always looks for f = 0).
//   - [t0, t1] is already clipped to the bounding box by the caller.
//   - Returns the smallest positive root in [t0, t1], or INFTY if none found.
//==============================================================================

class ImplicitSolver {
public:
  ImplicitSolver(double atol = 1e-8, int max_iter = 1000000)
    : atol_(atol), max_iter_(max_iter)
  {}
  virtual ~ImplicitSolver() = default;

  //! Find the smallest root of func(r + t*u) = 0 in [t0, t1].
  virtual double solve(const Implicit& func, Position r, Direction u, double t0,
    double t1, double isovalue = 0.0) const = 0;

  // access atol outside the solver
  double get_atol() const { return atol_; }

protected:
  double atol_;
  int max_iter_;
};

//==============================================================================
// NaiveLipschitz
//
// Lipschitz-marching solver.  At each step, the Lipschitz constant L gives a
// guaranteed safe skip distance: if |f(t)| / L > remaining interval, no root
// is possible and we advance.  When a sign change is detected, bisection
// finds the root to within atol.
//
// The Lipschitz constant is computed once over [t0, t1] at the start and
// reused throughout the march.  This is the "naive" part — a more
// sophisticated solver would recompute L on each sub-interval.
//
// Parameters:
//   atol     Geometric tolerance for root detection (metres). A value of
//            1e-8 (roughly one atomic radius) is a sensible default.
//   max_iter Safety cap on the number of marching steps.  If reached, the
//            solver returns INFTY and the surface is treated as not hit.
//==============================================================================

class NaiveLipschitz : public ImplicitSolver {
public:
  using ImplicitSolver::ImplicitSolver;
  double solve(const Implicit& function, Position r, Direction u, double t0,
    double t1, double isovalue) const override;

private:
  //! Bisect on [ta, tb] given that f(ta) and f(tb) have opposite signs.
  //! Returns the root to within atol_.
  double bisect(const Implicit& func, Position r, Direction u, double ta,
    double tb, double fa, double fb, double isovalue) const;
};

} // namespace openmc
#endif // OPENMC_IMPLICIT_SOLVER_H