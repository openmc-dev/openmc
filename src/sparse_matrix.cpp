//! \file sparse_matrix.cpp
//! \brief Implementation of CSCPattern and CSCMatrix

#include "openmc/sparse_matrix.h"

#include <algorithm>
#include <stdexcept>

#include <fmt/core.h>

namespace openmc {

//==============================================================================
// CSCPattern implementation
//==============================================================================

CSCPattern::CSCPattern(int n, vector<int> indptr, vector<int> indices)
  : n_(n), indptr_(std::move(indptr)), indices_(std::move(indices))
{
  if (static_cast<int>(indptr_.size()) != n_ + 1) {
    throw std::invalid_argument {fmt::format(
      "CSCPattern: indptr size ({}) != n + 1 ({})", indptr_.size(), n_ + 1)};
  }
  if (indptr_[0] != 0) {
    throw std::invalid_argument {
      fmt::format("CSCPattern: indptr[0] ({}) != 0", indptr_[0])};
  }
  if (indptr_[n_] != static_cast<int>(indices_.size())) {
    throw std::invalid_argument {
      fmt::format("CSCPattern: indptr[n] ({}) != indices size ({})",
        indptr_[n_], indices_.size())};
  }

  for (int j = 0; j < n_; ++j) {
    for (int k = indptr_[j]; k < indptr_[j + 1]; ++k) {
      if (indices_[k] < 0 || indices_[k] >= n_) {
        throw std::invalid_argument {fmt::format(
          "CSCPattern: row index {} out of bounds [0, {}) in column {}",
          indices_[k], n_, j)};
      }
      if (k > indptr_[j] && indices_[k - 1] >= indices_[k]) {
        throw std::invalid_argument {
          fmt::format("CSCPattern: row indices not sorted in column {} "
                      "(indices[{}]={} >= indices[{}]={})",
            j, k - 1, indices_[k - 1], k, indices_[k])};
      }
    }
  }
}

//==============================================================================
// CSCMatrix implementation
//==============================================================================

CSCMatrix::CSCMatrix(
  int n, vector<int> indptr, vector<int> indices, vector<double> data)
  : pattern_(n, std::move(indptr), std::move(indices)), data_(std::move(data))
{
  if (static_cast<int>(data_.size()) != pattern_.nnz()) {
    throw std::invalid_argument {
      fmt::format("CSCMatrix: data size ({}) != pattern nnz ({})",
        data_.size(), pattern_.nnz())};
  }
}

bool CSCPattern::operator==(const CSCPattern& other) const
{
  return n_ == other.n_ && indptr_ == other.indptr_ &&
         indices_ == other.indices_;
}

CSCPattern CSCPattern::with_diagonal() const
{
  // First pass: count entries per column, noting missing diagonals
  int extra = 0;
  for (int col = 0; col < n_; ++col) {
    bool has_diag = false;
    for (int idx = indptr_[col]; idx < indptr_[col + 1]; ++idx) {
      if (indices_[idx] == col) {
        has_diag = true;
        break;
      }
    }
    if (!has_diag)
      ++extra;
  }

  if (extra == 0)
    return *this;

  // Build new CSC directly, inserting diagonal entries in sorted position
  int new_nnz = nnz() + extra;
  vector<int> new_indptr(n_ + 1);
  vector<int> new_indices(new_nnz);

  int dst = 0;
  for (int col = 0; col < n_; ++col) {
    new_indptr[col] = dst;
    bool has_diag = false;
    bool diag_inserted = false;
    for (int idx = indptr_[col]; idx < indptr_[col + 1]; ++idx) {
      int row = indices_[idx];
      if (!has_diag && !diag_inserted && row > col) {
        new_indices[dst++] = col;
        diag_inserted = true;
      }
      if (row == col)
        has_diag = true;
      new_indices[dst++] = row;
    }
    if (!has_diag && !diag_inserted) {
      new_indices[dst++] = col;
    }
  }
  new_indptr[n_] = dst;

  return CSCPattern(n_, std::move(new_indptr), std::move(new_indices));
}

//==============================================================================
// LU factorization helpers
//==============================================================================

SymbolicLUFactorization symbolic_factorize(CSCPattern pattern)
{
  int n = pattern.n();
  const auto& indptr = pattern.indptr();
  const auto& indices = pattern.indices();

  // Build L and U fill patterns column by column.
  vector<vector<int>> l_cols(n);
  vector<bool> marked(n, false);
  vector<int> u_work, l_work;

  // Temporary storage for U column patterns before CSC assembly.
  vector<vector<int>> u_cols(n);

  for (int j = 0; j < n; ++j) {
    u_work.clear();
    l_work.clear();

    // Scatter: mark off-diagonal structural entries in column j.
    for (int p = indptr[j]; p < indptr[j + 1]; ++p) {
      int i = indices[p];
      if (i == j)
        continue;
      if (!marked[i]) {
        marked[i] = true;
        if (i < j) {
          u_work.push_back(i);
        } else {
          l_work.push_back(i);
        }
      }
    }

    // Propagate fill through previously discovered L columns.
    for (size_t idx = 0; idx < u_work.size(); ++idx) {
      int k = u_work[idx];
      for (int row : l_cols[k]) {
        if (row == j)
          continue;
        if (!marked[row]) {
          marked[row] = true;
          if (row < j) {
            u_work.push_back(row);
          } else {
            l_work.push_back(row);
          }
        }
      }
    }

    std::sort(u_work.begin(), u_work.end());
    std::sort(l_work.begin(), l_work.end());
    u_cols[j] = u_work;
    l_cols[j] = l_work;

    for (int k : u_work)
      marked[k] = false;
    for (int i : l_work)
      marked[i] = false;
  }

  vector<int> l_indptr(n + 1);
  vector<int> l_rowidx;
  for (int j = 0; j < n; ++j) {
    l_indptr[j] = static_cast<int>(l_rowidx.size());
    for (int r : l_cols[j])
      l_rowidx.push_back(r);
  }
  l_indptr[n] = static_cast<int>(l_rowidx.size());

  // Store the diagonal as the last U entry in each column.
  vector<int> u_indptr(n + 1);
  vector<int> u_rowidx;
  for (int j = 0; j < n; ++j) {
    u_indptr[j] = static_cast<int>(u_rowidx.size());
    for (int r : u_cols[j])
      u_rowidx.push_back(r);
    u_rowidx.push_back(j);
  }
  u_indptr[n] = static_cast<int>(u_rowidx.size());

  return {std::move(pattern),
    CSCPattern(n, std::move(l_indptr), std::move(l_rowidx)),
    CSCPattern(n, std::move(u_indptr), std::move(u_rowidx))};
}

} // namespace openmc
