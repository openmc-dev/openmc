//! \file bateman_solvers.h
//! \brief Solvers for the Bateman depletion equations

#ifndef OPENMC_BATEMAN_SOLVERS_H
#define OPENMC_BATEMAN_SOLVERS_H

#include <complex>

#include "openmc/sparse_matrix.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Numeric LU factorization values for a CRAM shifted linear system
//!
//! Stores the complex-valued L/U entries for (dt*A - theta*I) corresponding
//! to a previously computed SymbolicLUFactorization. The U diagonal is stored
//! both in the CSC values array and as precomputed reciprocals for faster back
//! substitution.
//==============================================================================

struct NumericLUFactorization {
  vector<std::complex<double>> l_data;
  vector<std::complex<double>> u_data;
  vector<std::complex<double>> u_diag_inv;
};

//! Numeric left-looking LU factorization of the CRAM shifted operator
//! (dt*A - theta*I) without explicitly forming it. Uses no pivoting; this is
//! safe because A is Metzler (non-negative off-diagonal) and |Im(theta)| >=
//! 1.194 for every CRAM pole, so every pivot satisfies |u_jj| >= |Im(theta)|.
void numeric_factorize_cram(const CSCMatrix& A, double dt,
  std::complex<double> theta, const SymbolicLUFactorization& symbolic,
  NumericLUFactorization& numeric, vector<std::complex<double>>& work);

//! Solve LUx = b using a symbolic LU pattern and matching numeric values.
void triangular_solve_lu(const vector<double>& b,
  const SymbolicLUFactorization& symbolic,
  const NumericLUFactorization& numeric, vector<std::complex<double>>& x);

//==============================================================================
//! Abstract base class for Bateman equation solvers
//!
//! Solves dN/dt = A*N over an interval dt, given initial composition N(0).
//==============================================================================

class BatemanSolver {
public:
  virtual ~BatemanSolver() = default;

  //! Solve the Bateman equations
  //! \param A   Sparse transmutation matrix (n x n)
  //! \param n0  Initial atom densities [n]
  //! \param dt  Time interval [s]
  //! \param substeps Number of substeps to use within dt
  //! \return    Final atom densities [n]
  virtual vector<double> solve(const CSCMatrix& A, const vector<double>& n0,
    double dt, int substeps = 1) = 0;
};

//==============================================================================
//! IPF CRAM solver for general transmutation matrices
//!
//! Implements the Incomplete Partial Fraction form of the Chebyshev
//! Rational Approximation Method (CRAM), as described in:
//!   M. Pusa, "Higher-Order Chebyshev Rational Approximation Method and
//!   Application to Burnup Equations," Nucl. Sci. Eng., 182:3, 297-318 (2016).
//!
//! The numeric factorization uses left-looking column LU without
//! pivoting. Each call to solve() performs a symbolic factorization
//! to compute L/U sparsity patterns, then reuses those patterns for
//! all pole solves within that call. Pivoting is unnecessary
//! because the transmutation matrix is Metzler (non-negative off-diagonal)
//! and the CRAM poles have nonzero imaginary parts (|Im(theta)| >= 1.194).
//! For any real Metzler matrix R and complex shift theta with Im(theta) != 0,
//! unpivoted Gaussian elimination on (R - theta*I) produces pivots u_jj
//! satisfying |u_jj| >= |Im(theta)|, guaranteeing non-singular factorization.
//! Since pivoting is not needed, the L/U sparsity patterns are deterministic
//! and identical across all poles, allowing a single symbolic phase to serve
//! all 24 (CRAM48) or 8 (CRAM16) linear solves.
//==============================================================================

class IPFCramSolver : public BatemanSolver {
public:
  explicit IPFCramSolver(int order = 48);

  //! Solve using full LU factorization (general transmutation matrix).
  vector<double> solve(const CSCMatrix& A, const vector<double>& n0, double dt,
    int substeps = 1) override;

private:
  // --- CRAM coefficients (non-owning views of the static tables) ---
  int n_poles_;                       //!< Number of poles (k/2)
  const std::complex<double>* alpha_; //!< Residues [n_poles]
  const std::complex<double>* theta_; //!< Poles [n_poles]
  double alpha0_;                     //!< Limit at infinity
};

} // namespace openmc

#endif // OPENMC_BATEMAN_SOLVERS_H
