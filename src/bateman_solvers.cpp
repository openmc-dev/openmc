//! \file bateman_solvers.cpp
//! \brief Implementation of Bateman equation solvers

#include "openmc/bateman_solvers.h"

#include <algorithm> // for sort
#include <complex>

#include "openmc/capi.h"
#include "openmc/error.h"

namespace openmc {

namespace {

// Fast complex reciprocal: 1/(a+bi) = (a-bi) / (a² + b²)
// Avoids GCC's __divdc3 which includes NaN/Inf handling we don't need.
// CRAM poles guarantee the diagonal is always well-conditioned.
inline std::complex<double> fast_crecip(std::complex<double> z)
{
  double a = z.real();
  double b = z.imag();
  double denom = a * a + b * b;
  return {a / denom, -b / denom};
}

// Fast complex multiply: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
// Avoids GCC's __muldc3 which includes NaN recovery we don't need.
inline std::complex<double> fast_cmul(
  std::complex<double> x, std::complex<double> y)
{
  double a = x.real(), b = x.imag();
  double c = y.real(), d = y.imag();
  return {a * c - b * d, a * d + b * c};
}

} // anonymous namespace

//==============================================================================
// CRAM coefficient tables
//
// Coefficients from M. Pusa, "Higher-Order Chebyshev Rational Approximation
// Method and Application to Burnup Equations," Nucl. Sci. Eng., 182:3,
// 297-318 (2016). Values match openmc/deplete/cram.py exactly.
//==============================================================================

namespace {

// --- CRAM16 coefficients (8 poles) ---

const std::complex<double> cram16_theta[] = {
  {+3.509103608414918e+0, +8.436198985884374e+0},
  {+5.948152268951177e+0, +3.587457362018322e+0},
  {-5.264971343442647e+0, +1.622022147316793e+1},
  {+1.419375897185666e+0, +1.092536348449672e+1},
  {+6.416177699099435e+0, +1.194122393370139e+0},
  {+4.993174737717997e+0, +5.996881713603942e+0},
  {-1.413928462488886e+0, +1.349772569889275e+1},
  {-1.084391707869699e+1, +1.927744616718165e+1},
};

const std::complex<double> cram16_alpha[] = {
  {+5.464930576870210e+3, -3.797983575308356e+4},
  {+9.045112476907548e+1, -1.115537522430261e+3},
  {+2.344818070467641e+2, -4.228020157070496e+2},
  {+9.453304067358312e+1, -2.951294291446048e+2},
  {+7.283792954673409e+2, -1.205646080220011e+5},
  {+3.648229059594851e+1, -1.155509621409682e+2},
  {+2.547321630156819e+1, -2.639500283021502e+1},
  {+2.394538338734709e+1, -5.650522971778156e+0},
};

const double cram16_alpha0 = 2.124853710495224e-16;

// --- CRAM48 coefficients (24 poles) ---

const std::complex<double> cram48_theta[] = {
  {-4.465731934165702e+1, +6.233225190695437e+1},
  {-5.284616241568964e+0, +4.057499381311059e+1},
  {-8.867715667624458e+0, +4.325515754166724e+1},
  {+3.493013124279215e+0, +3.281615453173585e+1},
  {+1.564102508858634e+1, +1.558061616372237e+1},
  {+1.742097597385893e+1, +1.076629305714420e+1},
  {-2.834466755180654e+1, +5.492841024648724e+1},
  {+1.661569367939544e+1, +1.316994930024688e+1},
  {+8.011836167974721e+0, +2.780232111309410e+1},
  {-2.056267541998229e+0, +3.794824788914354e+1},
  {+1.449208170441839e+1, +1.799988210051809e+1},
  {+1.853807176907916e+1, +5.974332563100539e+0},
  {+9.932562704505182e+0, +2.532823409972962e+1},
  {-2.244223871767187e+1, +5.179633600312162e+1},
  {+8.590014121680897e-1, +3.536456194294350e+1},
  {-1.286192925744479e+1, +4.600304902833652e+1},
  {+1.164596909542055e+1, +2.287153304140217e+1},
  {+1.806076684783089e+1, +8.368200580099821e+0},
  {+5.870672154659249e+0, +3.029700159040121e+1},
  {-3.542938819659747e+1, +5.834381701800013e+1},
  {+1.901323489060250e+1, +1.194282058271408e+0},
  {+1.885508331552577e+1, +3.583428564427879e+0},
  {-1.734689708174982e+1, +4.883941101108207e+1},
  {+1.316284237125190e+1, +2.042951874827759e+1},
};

const std::complex<double> cram48_alpha[] = {
  {+6.387380733878774e+2, -6.743912502859256e+2},
  {+1.909896179065730e+2, -3.973203432721332e+2},
  {+4.236195226571914e+2, -2.041233768918671e+3},
  {+4.645770595258726e+2, -1.652917287299683e+3},
  {+7.765163276752433e+2, -1.783617639907328e+4},
  {+1.907115136768522e+3, -5.887068595142284e+4},
  {+2.909892685603256e+3, -9.953255345514560e+3},
  {+1.944772206620450e+2, -1.427131226068449e+3},
  {+1.382799786972332e+5, -3.256885197214938e+6},
  {+5.628442079602433e+3, -2.924284515884309e+4},
  {+2.151681283794220e+2, -1.121774011188224e+3},
  {+1.324720240514420e+3, -6.370088443140973e+4},
  {+1.617548476343347e+4, -1.008798413156542e+6},
  {+1.112729040439685e+2, -8.837109731680418e+1},
  {+1.074624783191125e+2, -1.457246116408180e+2},
  {+8.835727765158191e+1, -6.388286188419360e+1},
  {+9.354078136054179e+1, -2.195424319460237e+2},
  {+9.418142823531573e+1, -6.719055740098035e+2},
  {+1.040012390717851e+2, -1.693747595553868e+2},
  {+6.861882624343235e+1, -1.177598523430493e+1},
  {+8.766654491283722e+1, -4.596464999363902e+3},
  {+1.056007619389650e+2, -1.738294585524067e+3},
  {+7.738987569039419e+1, -4.311715386228984e+1},
  {+1.041366366475571e+2, -2.777743732451969e+2},
};

const double cram48_alpha0 = 2.258038182743983e-47;

} // anonymous namespace

//==============================================================================
// IPFCramSolver implementation
//==============================================================================

IPFCramSolver::IPFCramSolver(CramOrder order)
{
  if (order == CramOrder::cram16) {
    n_poles_ = 8;
    alpha_.assign(cram16_alpha, cram16_alpha + n_poles_);
    theta_.assign(cram16_theta, cram16_theta + n_poles_);
    alpha0_ = cram16_alpha0;
  } else {
    n_poles_ = 24;
    alpha_.assign(cram48_alpha, cram48_alpha + n_poles_);
    theta_.assign(cram48_theta, cram48_theta + n_poles_);
    alpha0_ = cram48_alpha0;
  }
}

vector<double> IPFCramSolver::solve(
  const CSCMatrix& A, const vector<double>& n0, double dt)
{
  int n = A.n();

  // Symbolic factorization: compute L/U sparsity patterns for this matrix.
  // Reused across all pole solves below.
  CSCPattern pattern = A.pattern().with_diagonal();
  symbolic_factorize(pattern);

  // IPF CRAM iteration:
  //   y_0 = n0
  //   y_{k+1} = y_k + 2*Re(alpha_k * (A*dt - theta_k*I)^{-1} * y_k)
  //   result = alpha0 * y_final
  vector<double> y(n0.begin(), n0.end());

  for (int p = 0; p < n_poles_; ++p) {
    numeric_factorize(A, pattern, dt, theta_[p]);
    triangular_solve(y, x_);

    // y += 2 * Re(alpha_p * x_p)
    for (int i = 0; i < n; ++i) {
      auto ax = fast_cmul(alpha_[p], x_[i]);
      y[i] += 2.0 * ax.real();
    }
  }

  // Final scaling
  for (int i = 0; i < n; ++i) {
    y[i] *= alpha0_;
  }

  return y;
}

//==============================================================================
// Symbolic factorization
//
// Computes the exact L/U sparsity patterns for left-looking column LU
// factorization without pivoting. Pivoting is unnecessary because the
// transmutation matrix A is Metzler (non-negative off-diagonal entries)
// and each CRAM pole theta has nonzero imaginary part. For M = A*dt - theta*I
// with A Metzler, unpivoted Gaussian elimination produces pivots u_jj
// satisfying |u_jj| >= |Im(theta)| >= 1.194. This guarantees non-singular
// factorization, and since no row swaps occur the L/U patterns are
// deterministic and identical across all poles.
//
// Algorithm: symbolic left-looking factorization with worklist-based fill
// propagation. For each column j, start with the structural nonzeros of
// the input pattern, then propagate fill through previously computed L
// column patterns. Any row k < j that becomes nonzero (a U entry) triggers
// examination of L[:,k]'s rows, which may create additional fill.
//==============================================================================

void IPFCramSolver::symbolic_factorize(const CSCPattern& pattern)
{
  int n = pattern.n();
  const auto& indptr = pattern.indptr();
  const auto& indices = pattern.indices();

  // Build L and U fill patterns column by column
  vector<vector<int>> l_cols(n);
  vector<bool> marked(n, false);
  vector<int> u_work, l_work;

  // Temporary storage for U column patterns (needed for CSC construction)
  vector<vector<int>> u_cols(n);

  for (int j = 0; j < n; ++j) {
    u_work.clear();
    l_work.clear();

    // Scatter: mark off-diagonal nonzero rows of column j
    for (int p = indptr[j]; p < indptr[j + 1]; ++p) {
      int i = indices[p];
      if (i == j)
        continue;
      if (!marked[i]) {
        marked[i] = true;
        if (i < j)
          u_work.push_back(i);
        else
          l_work.push_back(i);
      }
    }

    // Propagate fill through L columns of above-diagonal entries.
    // u_work grows as new above-diagonal rows are discovered via fill.
    for (size_t idx = 0; idx < u_work.size(); ++idx) {
      int k = u_work[idx];
      for (int row : l_cols[k]) {
        if (row == j)
          continue; // diagonal, skip
        if (!marked[row]) {
          marked[row] = true;
          if (row < j)
            u_work.push_back(row);
          else
            l_work.push_back(row);
        }
      }
    }

    // Sort row indices and store
    std::sort(u_work.begin(), u_work.end());
    std::sort(l_work.begin(), l_work.end());
    u_cols[j] = u_work;
    l_cols[j] = l_work;

    // Clear marked flags
    for (int k : u_work)
      marked[k] = false;
    for (int i : l_work)
      marked[i] = false;
  }

  // Build L CSC-style index arrays
  l_indptr_.resize(n + 1);
  l_rowidx_.clear();
  for (int j = 0; j < n; ++j) {
    l_indptr_[j] = static_cast<int>(l_rowidx_.size());
    for (int r : l_cols[j])
      l_rowidx_.push_back(r);
  }
  l_indptr_[n] = static_cast<int>(l_rowidx_.size());

  // Build U CSC-style index arrays (diagonal stored as last entry per column)
  u_indptr_.resize(n + 1);
  u_rowidx_.clear();
  for (int j = 0; j < n; ++j) {
    u_indptr_[j] = static_cast<int>(u_rowidx_.size());
    for (int r : u_cols[j])
      u_rowidx_.push_back(r);
    u_rowidx_.push_back(j); // diagonal last
  }
  u_indptr_[n] = static_cast<int>(u_rowidx_.size());

  // Allocate numeric workspace
  l_data_.resize(l_rowidx_.size());
  u_data_.resize(u_rowidx_.size());
  u_diag_.resize(n);
  work_.resize(n);
  x_.resize(n);
}

//==============================================================================
// Numeric factorization
//
// Left-looking column LU factorization without pivoting.
// Forms the shifted matrix M = A*dt - theta*I on-the-fly.
//
// For each column j, the above-diagonal U rows (stored in ascending order
// in u_rowidx_) serve as the left-looking schedule: each k < j with
// U[k,j] != 0 triggers a rank-1 update w -= L[:,k] * U[k,j]. Processing
// in ascending order naturally performs the forward substitution that
// resolves fill-in dependencies between earlier columns.
//==============================================================================

void IPFCramSolver::numeric_factorize(const CSCMatrix& A,
  const CSCPattern& pattern, double dt, std::complex<double> theta)
{
  int n = A.n();
  const auto& a_indptr = A.indptr();
  const auto& a_indices = A.indices();
  const auto& a_data = A.data();

  const auto& sp_indptr = pattern.indptr();
  const auto& sp_indices = pattern.indices();

  for (int j = 0; j < n; ++j) {

    // --- Step 1: Scatter M[:,j] = A[:,j]*dt - theta*I[:,j] ---
    {
      int a_pos = a_indptr[j];
      int a_end = a_indptr[j + 1];

      for (int sp_pos = sp_indptr[j]; sp_pos < sp_indptr[j + 1]; ++sp_pos) {
        int row = sp_indices[sp_pos];
        std::complex<double> val(0.0, 0.0);

        if (a_pos < a_end && a_indices[a_pos] == row) {
          val = dt * a_data[a_pos];
          ++a_pos;
        }

        if (row == j) {
          val -= theta;
        }

        work_[row] = val;
      }
    }

    // --- Step 2: Left-looking updates ---
    // Process U[:,j] rows in ascending order (forward substitution).
    // Each k < j with U[k,j] != 0 subtracts L[:,k] * U[k,j].
    for (int up = u_indptr_[j]; up < u_indptr_[j + 1] - 1; ++up) {
      int k = u_rowidx_[up];
      std::complex<double> ukj = work_[k];
      u_data_[up] = ukj;

      for (int lp = l_indptr_[k]; lp < l_indptr_[k + 1]; ++lp) {
        work_[l_rowidx_[lp]] -= fast_cmul(l_data_[lp], ukj);
      }
    }

    // --- Step 3: Extract diagonal and L column ---
    std::complex<double> inv_ujj = fast_crecip(work_[j]);
    u_diag_[j] = inv_ujj;
    u_data_[u_indptr_[j + 1] - 1] = work_[j];
    for (int lp = l_indptr_[j]; lp < l_indptr_[j + 1]; ++lp) {
      l_data_[lp] = fast_cmul(work_[l_rowidx_[lp]], inv_ujj);
    }

    // --- Step 4: Clear workspace ---
    // Clear all positions in the predicted fill pattern (U + L + diagonal)
    for (int up = u_indptr_[j]; up < u_indptr_[j + 1]; ++up) {
      work_[u_rowidx_[up]] = {0.0, 0.0};
    }
    for (int lp = l_indptr_[j]; lp < l_indptr_[j + 1]; ++lp) {
      work_[l_rowidx_[lp]] = {0.0, 0.0};
    }
  }
}

//==============================================================================
// Triangular solve
//
// Solve LUx = b without permutation (no pivoting means identity perm).
// Forward substitution solves Lz = b, back substitution solves Ux = z.
//==============================================================================

void IPFCramSolver::triangular_solve(
  const vector<double>& b, vector<std::complex<double>>& x) const
{
  int n = static_cast<int>(u_diag_.size());

  // Copy real RHS into complex vector
  for (int j = 0; j < n; ++j) {
    x[j] = std::complex<double>(b[j], 0.0);
  }

  // Forward substitution: Lz = b (L is unit lower triangular)
  for (int j = 0; j < n; ++j) {
    for (int lp = l_indptr_[j]; lp < l_indptr_[j + 1]; ++lp) {
      x[l_rowidx_[lp]] -= fast_cmul(l_data_[lp], x[j]);
    }
  }

  // Back substitution: Ux = z
  // u_diag_ stores reciprocals (1/U[j,j]) to avoid complex division
  for (int j = n - 1; j >= 0; --j) {
    x[j] = fast_cmul(x[j], u_diag_[j]);

    for (int up = u_indptr_[j]; up < u_indptr_[j + 1] - 1; ++up) {
      x[u_rowidx_[up]] -= fast_cmul(u_data_[up], x[j]);
    }
  }
}

} // namespace openmc

//==============================================================================
// C API
//==============================================================================

using namespace openmc;

extern "C" int openmc_cram_solve(int n, const int* indptr, const int* indices,
  const double* data, const double* n0, double dt, int order, double* result)
{
  try {
    if (order != 16 && order != 48) {
      set_errmsg(fmt::format("CRAM order must be 16 or 48, got {}", order));
      return OPENMC_E_INVALID_ARGUMENT;
    }

    auto cram_order = (order == 16) ? CramOrder::cram16 : CramOrder::cram48;

    int nnz = indptr[n];
    CSCMatrix A(n, vector<int>(indptr, indptr + n + 1),
      vector<int>(indices, indices + nnz), vector<double>(data, data + nnz));
    vector<double> n0_vec(n0, n0 + n);

    IPFCramSolver solver(cram_order);
    vector<double> y = solver.solve(A, n0_vec, dt);
    std::copy(y.begin(), y.end(), result);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}
