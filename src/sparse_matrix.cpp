//! \file sparse_matrix.cpp
//! \brief Implementation of CSCPattern and CSCMatrix

#include "openmc/sparse_matrix.h"

#include "openmc/error.h"

namespace openmc {

//==============================================================================
// CSCPattern implementation
//==============================================================================

CSCPattern::CSCPattern(int n, vector<int> indptr, vector<int> indices)
  : n_(n), indptr_(std::move(indptr)), indices_(std::move(indices))
{
  if (static_cast<int>(indptr_.size()) != n_ + 1) {
    fatal_error(fmt::format(
      "CSCPattern: indptr size ({}) != n + 1 ({})", indptr_.size(), n_ + 1));
  }
  if (indptr_[0] != 0) {
    fatal_error(fmt::format("CSCPattern: indptr[0] ({}) != 0", indptr_[0]));
  }
  if (indptr_[n_] != static_cast<int>(indices_.size())) {
    fatal_error(fmt::format(
      "CSCPattern: indptr[n] ({}) != indices size ({})",
      indptr_[n_], indices_.size()));
  }

  for (int j = 0; j < n_; ++j) {
    for (int k = indptr_[j]; k < indptr_[j + 1]; ++k) {
      if (indices_[k] < 0 || indices_[k] >= n_) {
        fatal_error(fmt::format(
          "CSCPattern: row index {} out of bounds [0, {}) in column {}",
          indices_[k], n_, j));
      }
      if (k > indptr_[j] && indices_[k - 1] >= indices_[k]) {
        fatal_error(fmt::format(
          "CSCPattern: row indices not sorted in column {} "
          "(indices[{}]={} >= indices[{}]={})",
          j, k - 1, indices_[k - 1], k, indices_[k]));
      }
    }
  }
}

//==============================================================================
// CSCMatrix implementation
//==============================================================================

CSCMatrix::CSCMatrix(
  int n, vector<int> indptr, vector<int> indices, vector<double> data)
  : pattern_(n, std::move(indptr), std::move(indices)),
    data_(std::move(data))
{
  if (static_cast<int>(data_.size()) != pattern_.nnz()) {
    fatal_error(fmt::format("CSCMatrix: data size ({}) != pattern nnz ({})",
      data_.size(), pattern_.nnz()));
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

} // namespace openmc
