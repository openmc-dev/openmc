//! \file sparse_matrix.h
//! \brief Compressed Sparse Column (CSC) sparsity pattern and matrix classes

#ifndef OPENMC_SPARSE_MATRIX_H
#define OPENMC_SPARSE_MATRIX_H

#include <utility> // for move

#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! CSC sparsity pattern (structure only, no values)
//!
//! Stores the compressed column structure of a square sparse matrix:
//! column pointers and row indices. Rows within each column are sorted in
//! ascending order.
//==============================================================================

class CSCPattern {
public:
  // Constructors
  CSCPattern() = default;
  CSCPattern(int n, vector<int> indptr, vector<int> indices)
    : n_(n), indptr_(std::move(indptr)), indices_(std::move(indices))
  {}

  // Accessors
  int n() const { return n_; }
  int nnz() const { return static_cast<int>(indices_.size()); }
  const vector<int>& indptr() const { return indptr_; }
  const vector<int>& indices() const { return indices_; }

  //! Return a new pattern with all diagonal entries forced present.
  //! Existing entries (including any diagonals already present) are preserved.
  CSCPattern with_diagonal() const;

  //! Structural equality check. Two patterns are equal iff they have the same
  //! dimension, identical column pointers, and identical row indices.
  bool operator==(const CSCPattern& other) const;
  bool operator!=(const CSCPattern& other) const { return !(*this == other); }

private:
  int n_ {0};           //!< Matrix dimension
  vector<int> indptr_;  //!< Column pointers [n+1]
  vector<int> indices_; //!< Row indices [nnz], sorted within each column
};

//==============================================================================
//! CSC sparse matrix with real (double) values
//!
//! Associates a double-precision value with each structural nonzero defined
//! by the underlying CSCPattern.
//==============================================================================

class CSCMatrix {
public:
  // Constructors
  CSCMatrix() = default;
  CSCMatrix(CSCPattern pattern, vector<double> data)
    : pattern_(std::move(pattern)), data_(std::move(data))
  {}

  // Accessors
  int n() const { return pattern_.n(); }
  int nnz() const { return pattern_.nnz(); }
  const CSCPattern& pattern() const { return pattern_; }
  const vector<int>& indptr() const { return pattern_.indptr(); }
  const vector<int>& indices() const { return pattern_.indices(); }
  const vector<double>& data() const { return data_; }

private:
  CSCPattern pattern_;  //!< Structural pattern
  vector<double> data_; //!< Values [nnz]
};

} // namespace openmc

#endif // OPENMC_SPARSE_MATRIX_H
