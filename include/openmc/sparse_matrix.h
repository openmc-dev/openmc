//! \file sparse_matrix.h
//! \brief Compressed Sparse Column (CSC) sparsity pattern and matrix classes

#ifndef OPENMC_SPARSE_MATRIX_H
#define OPENMC_SPARSE_MATRIX_H

#include <algorithm> // for sort, adjacent_find
#include <utility>   // for move

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

  //! Construct from coordinate (COO) triplets.
  //! Duplicate (row, col) pairs are allowed — only the structural pattern
  //! is retained (duplicates are collapsed to a single entry).
  //! \param n Matrix dimension (square n x n)
  //! \param rows Row indices
  //! \param cols Column indices
  static CSCPattern from_triplets(
    int n, const vector<int>& rows, const vector<int>& cols);

  // Accessors
  int n() const { return n_; }
  int nnz() const { return static_cast<int>(indices_.size()); }
  const vector<int>& indptr() const { return indptr_; }
  const vector<int>& indices() const { return indices_; }

  //! Return a new pattern with rows and columns permuted.
  //! \param perm Permutation vector: new_index -> old_index
  CSCPattern permute(const vector<int>& perm) const;

  //! Return a new pattern with all diagonal entries forced present.
  //! Existing entries (including any diagonals already present) are preserved.
  CSCPattern with_diagonal() const;

  //! Structural equality check. Two patterns are equal iff they have the same
  //! dimension, identical column pointers, and identical row indices.
  bool operator==(const CSCPattern& other) const;
  bool operator!=(const CSCPattern& other) const { return !(*this == other); }

  //! Compute structural reachability (transitive closure) under a topological
  //! permutation. For each column j in permuted space, computes the set of row
  //! indices reachable via the directed graph defined by off-diagonal entries.
  //! The result is invariant for a given topology and can be reused across
  //! time steps. Output is in flat CSC-like format.
  //!
  //! \param perm           Topological permutation: perm[new_idx] = old_idx
  //! \param reach_indptr   Output column pointers [n+1]
  //! \param reach_indices  Output row indices [total_reach]
  void reachability(const vector<int>& perm, vector<int>& reach_indptr,
    vector<int>& reach_indices) const;

  //! Compute a topological sort of the directed graph defined by the
  //! off-diagonal entries of this pattern. The graph has edge col -> row
  //! for each off-diagonal entry (row, col). Returns a permutation vector
  //! perm where perm[new_idx] = old_idx such that the permuted matrix is
  //! lower-triangular (rows >= col for all entries).
  //!
  //! Uses Kahn's algorithm (BFS peeling of zero-in-degree nodes). Nodes
  //! with equal in-degree are processed in ascending index order for
  //! determinism.
  //!
  //! \return Topological permutation vector of length n
  //! \throws std::runtime_error if the graph contains a cycle
  vector<int> topological_sort() const;

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

  //! Construct an empty (zero) n x n matrix.
  explicit CSCMatrix(int n)
    : pattern_(n, vector<int>(n + 1, 0), {}), data_()
  {}

  //! Construct from coordinate (COO) triplets.
  //! Duplicate (row, col) pairs are summed.
  //! \param n Matrix dimension (square n x n)
  //! \param rows Row indices
  //! \param cols Column indices
  //! \param vals Values
  static CSCMatrix from_triplets(int n, const vector<int>& rows,
    const vector<int>& cols, const vector<double>& vals);

  // Accessors
  int n() const { return pattern_.n(); }
  int nnz() const { return pattern_.nnz(); }
  const CSCPattern& pattern() const { return pattern_; }
  const vector<int>& indptr() const { return pattern_.indptr(); }
  const vector<int>& indices() const { return pattern_.indices(); }
  const vector<double>& data() const { return data_; }

  //! Return a new matrix with rows and columns permuted.
  //! \param perm Permutation vector: new_index -> old_index
  CSCMatrix permute(const vector<int>& perm) const;

  //! Element-wise addition of two CSC matrices with the same dimension.
  //! The sparsity patterns may differ; the result has the union of both.
  CSCMatrix operator+(const CSCMatrix& other) const;

  //! In-place element-wise addition. Same semantics as operator+ but
  //! avoids allocating a new matrix when the sparsity pattern is a
  //! superset of \p other's pattern. Falls back to operator+ otherwise.
  CSCMatrix& operator+=(const CSCMatrix& other);

  //! Scalar multiplication (returns a new matrix).
  friend CSCMatrix operator*(double scalar, const CSCMatrix& mat);

  //! Scale all values in-place by \p scalar.
  void scale(double scalar);

private:
  CSCPattern pattern_;  //!< Structural pattern
  vector<double> data_; //!< Values [nnz]
};

} // namespace openmc

#endif // OPENMC_SPARSE_MATRIX_H
