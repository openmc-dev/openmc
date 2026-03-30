//! \file sparse_matrix.cpp
//! \brief Implementation of CSCPattern and CSCMatrix

#include "openmc/sparse_matrix.h"

#include <algorithm>  // for sort, fill
#include <functional> // for plus
#include <numeric>    // for iota
#include <utility>    // for pair

#include <fmt/core.h>

#include "openmc/error.h"

namespace openmc {

//==============================================================================
// CSCPattern implementation
//==============================================================================

CSCPattern CSCPattern::from_triplets(
  int n, const vector<int>& rows, const vector<int>& cols)
{
  int nt = rows.size();

  // Sort triplets by (col, row)
  vector<int> order(nt);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return cols[a] < cols[b] || (cols[a] == cols[b] && rows[a] < rows[b]);
  });

  // Build CSC arrays, collapsing duplicate (row, col) pairs
  vector<int> indptr(n + 1, 0);
  vector<int> indices;
  indices.reserve(nt);

  int prev_col = -1;
  int prev_row = -1;
  for (int k = 0; k < nt; ++k) {
    int i = rows[order[k]];
    int j = cols[order[k]];

    // Skip duplicate entries
    if (j == prev_col && i == prev_row)
      continue;

    indices.push_back(i);

    // Fill column pointers for any skipped columns
    for (int c = prev_col + 1; c <= j; ++c) {
      indptr[c] = static_cast<int>(indices.size()) - 1;
    }
    prev_col = j;
    prev_row = i;
  }
  // Fill remaining column pointers
  for (int c = prev_col + 1; c <= n; ++c) {
    indptr[c] = static_cast<int>(indices.size());
  }

  return CSCPattern(n, std::move(indptr), std::move(indices));
}

CSCPattern CSCPattern::permute(const vector<int>& perm) const
{
  // perm[new_index] = old_index
  // Build inverse: inv_perm[old_index] = new_index
  int n = n_;
  vector<int> inv_perm(n);
  for (int i = 0; i < n; ++i) {
    inv_perm[perm[i]] = i;
  }

  // Collect permuted triplets
  vector<int> new_rows, new_cols;
  for (int old_col = 0; old_col < n; ++old_col) {
    int new_col = inv_perm[old_col];
    for (int idx = indptr_[old_col]; idx < indptr_[old_col + 1]; ++idx) {
      int old_row = indices_[idx];
      int new_row = inv_perm[old_row];
      new_rows.push_back(new_row);
      new_cols.push_back(new_col);
    }
  }

  return CSCPattern::from_triplets(n, new_rows, new_cols);
}

bool CSCPattern::operator==(const CSCPattern& other) const
{
  return n_ == other.n_ && indptr_ == other.indptr_ &&
         indices_ == other.indices_;
}

void CSCPattern::reachability(const vector<int>& perm,
  vector<int>& reach_indptr, vector<int>& reach_indices) const
{
  int n = n_;

  // Build inverse permutation
  vector<int> inv_perm(n);
  for (int i = 0; i < n; ++i) {
    inv_perm[perm[i]] = i;
  }

  // Build lower-triangular column structure (off-diagonal only)
  vector<int> lt_indptr(n + 1, 0);
  for (int old_col = 0; old_col < n; ++old_col) {
    int new_col = inv_perm[old_col];
    for (int p = indptr_[old_col]; p < indptr_[old_col + 1]; ++p) {
      int new_row = inv_perm[indices_[p]];
      if (new_row != new_col) {
        ++lt_indptr[new_col + 1];
      }
    }
  }
  for (int j = 0; j < n; ++j) {
    lt_indptr[j + 1] += lt_indptr[j];
  }

  int lt_nnz = lt_indptr[n];
  vector<int> lt_rowidx(lt_nnz);
  vector<int> col_pos(n, 0);
  for (int old_col = 0; old_col < n; ++old_col) {
    int new_col = inv_perm[old_col];
    for (int p = indptr_[old_col]; p < indptr_[old_col + 1]; ++p) {
      int new_row = inv_perm[indices_[p]];
      if (new_row != new_col) {
        lt_rowidx[lt_indptr[new_col] + col_pos[new_col]++] = new_row;
      }
    }
  }

  // Sort row indices within each column
  for (int j = 0; j < n; ++j) {
    std::sort(lt_rowidx.begin() + lt_indptr[j],
      lt_rowidx.begin() + lt_indptr[j + 1]);
  }

  // Compute reach via memoized transitive closure (leaves first).
  // reach[j] = sorted indices of all nodes reachable from j (excluding j).
  vector<vector<int>> reach(n);
  vector<int> merge_buf;

  for (int j = n - 1; j >= 0; --j) {
    int n_children = lt_indptr[j + 1] - lt_indptr[j];
    if (n_children == 0)
      continue;

    if (n_children == 1) {
      int c = lt_rowidx[lt_indptr[j]];
      auto& rc = reach[c];
      reach[j].resize(1 + rc.size());
      auto it = std::lower_bound(rc.begin(), rc.end(), c);
      size_t pos = it - rc.begin();
      std::copy(rc.begin(), it, reach[j].begin());
      reach[j][pos] = c;
      std::copy(it, rc.end(), reach[j].begin() + pos + 1);
    } else {
      int c0 = lt_rowidx[lt_indptr[j]];
      auto& rc0 = reach[c0];
      merge_buf.clear();
      merge_buf.reserve(rc0.size() + 1);
      auto it0 = std::lower_bound(rc0.begin(), rc0.end(), c0);
      merge_buf.insert(merge_buf.end(), rc0.begin(), it0);
      merge_buf.push_back(c0);
      merge_buf.insert(merge_buf.end(), it0, rc0.end());

      for (int lp = lt_indptr[j] + 1; lp < lt_indptr[j + 1]; ++lp) {
        int c = lt_rowidx[lp];
        auto& rc = reach[c];
        vector<int> child_set;
        child_set.reserve(rc.size() + 1);
        auto itc = std::lower_bound(rc.begin(), rc.end(), c);
        child_set.insert(child_set.end(), rc.begin(), itc);
        child_set.push_back(c);
        child_set.insert(child_set.end(), itc, rc.end());

        vector<int> merged;
        merged.reserve(merge_buf.size() + child_set.size());
        std::set_union(merge_buf.begin(), merge_buf.end(), child_set.begin(),
          child_set.end(), std::back_inserter(merged));
        merge_buf = std::move(merged);
      }

      reach[j] = std::move(merge_buf);
    }
  }

  // Flatten to CSC-like format
  reach_indptr.resize(n + 1);
  reach_indptr[0] = 0;
  for (int j = 0; j < n; ++j) {
    reach_indptr[j + 1] =
      reach_indptr[j] + static_cast<int>(reach[j].size());
  }
  int total = reach_indptr[n];
  reach_indices.resize(total);
  for (int j = 0; j < n; ++j) {
    std::copy(
      reach[j].begin(), reach[j].end(), reach_indices.begin() + reach_indptr[j]);
  }
}

vector<int> CSCPattern::topological_sort() const
{
  int n = n_;

  // Build in-degree count from off-diagonal entries.
  // Graph edge: col -> row for each off-diagonal (row, col) entry.
  vector<int> in_degree(n, 0);
  for (int col = 0; col < n; ++col) {
    for (int p = indptr_[col]; p < indptr_[col + 1]; ++p) {
      int row = indices_[p];
      if (row != col) {
        ++in_degree[row];
      }
    }
  }

  // Initialize queue with zero-in-degree nodes (ascending order)
  // Using a min-heap for deterministic ordering.
  vector<int> queue;
  for (int i = 0; i < n; ++i) {
    if (in_degree[i] == 0) {
      queue.push_back(i);
    }
  }
  // Make a min-heap (smallest index first)
  std::make_heap(queue.begin(), queue.end(), std::greater<int>());

  vector<int> perm;
  perm.reserve(n);

  while (!queue.empty()) {
    std::pop_heap(queue.begin(), queue.end(), std::greater<int>());
    int node = queue.back();
    queue.pop_back();
    perm.push_back(node);

    // "Remove" node: decrement in-degree of successors
    for (int p = indptr_[node]; p < indptr_[node + 1]; ++p) {
      int row = indices_[p];
      if (row != node) {
        if (--in_degree[row] == 0) {
          queue.push_back(row);
          std::push_heap(queue.begin(), queue.end(), std::greater<int>());
        }
      }
    }
  }

  if (static_cast<int>(perm.size()) != n) {
    fatal_error(fmt::format(
      "Topological sort failed: graph contains a cycle "
      "({} of {} nodes processed)", perm.size(), n));
  }

  return perm;
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

  if (extra == 0) {
    return CSCPattern(n_, vector<int>(indptr_), vector<int>(indices_));
  }

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
// CSCMatrix implementation
//==============================================================================

CSCMatrix CSCMatrix::from_triplets(int n, const vector<int>& rows,
  const vector<int>& cols, const vector<double>& vals)
{
  int nt = rows.size();

  // Sort triplets by (col, row)
  vector<int> order(nt);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return cols[a] < cols[b] || (cols[a] == cols[b] && rows[a] < rows[b]);
  });

  // Build CSC arrays, summing duplicate (row, col) pairs
  vector<int> indptr(n + 1, 0);
  vector<int> indices;
  vector<double> data;
  indices.reserve(nt);
  data.reserve(nt);

  // First pass: group duplicates and sum values into a temporary buffer.
  // We need the sums before emitting entries so that zero sums are dropped.
  struct Entry {
    int row, col;
    double val;
  };
  vector<Entry> entries;
  entries.reserve(nt);

  int prev_col = -1;
  int prev_row = -1;
  for (int k = 0; k < nt; ++k) {
    int i = rows[order[k]];
    int j = cols[order[k]];
    double v = vals[order[k]];

    if (j == prev_col && i == prev_row) {
      entries.back().val += v;
      continue;
    }
    entries.push_back({i, j, v});
    prev_col = j;
    prev_row = i;
  }

  // Second pass: emit only nonzero entries
  prev_col = -1;
  for (auto& e : entries) {
    if (e.val == 0.0)
      continue;

    indices.push_back(e.row);
    data.push_back(e.val);

    for (int c = prev_col + 1; c <= e.col; ++c) {
      indptr[c] = static_cast<int>(indices.size()) - 1;
    }
    prev_col = e.col;
  }
  // Fill remaining column pointers
  for (int c = prev_col + 1; c <= n; ++c) {
    indptr[c] = static_cast<int>(indices.size());
  }

  CSCPattern pattern(n, std::move(indptr), std::move(indices));
  return CSCMatrix(std::move(pattern), std::move(data));
}

CSCMatrix CSCMatrix::permute(const vector<int>& perm) const
{
  // perm[new_index] = old_index
  int n = pattern_.n();
  vector<int> inv_perm(n);
  for (int i = 0; i < n; ++i) {
    inv_perm[perm[i]] = i;
  }

  // Collect permuted triplets with values
  vector<int> new_rows, new_cols;
  vector<double> new_vals;
  const auto& indptr = pattern_.indptr();
  const auto& indices = pattern_.indices();

  for (int old_col = 0; old_col < n; ++old_col) {
    int new_col = inv_perm[old_col];
    for (int idx = indptr[old_col]; idx < indptr[old_col + 1]; ++idx) {
      int old_row = indices[idx];
      int new_row = inv_perm[old_row];
      new_rows.push_back(new_row);
      new_cols.push_back(new_col);
      new_vals.push_back(data_[idx]);
    }
  }

  return CSCMatrix::from_triplets(n, new_rows, new_cols, new_vals);
}

CSCMatrix CSCMatrix::operator+(const CSCMatrix& other) const
{
  int n = pattern_.n();
  if (other.n() != n) {
    fatal_error(fmt::format(
      "Cannot add CSC matrices with different dimensions ({} vs {})", n,
      other.n()));
  }

  const auto& a_indptr = indptr();
  const auto& a_indices = indices();
  const auto& b_indptr = other.indptr();
  const auto& b_indices = other.indices();

  // Merge sorted row indices column-by-column
  vector<int> new_indptr(n + 1);
  vector<int> new_indices;
  vector<double> new_data;
  new_indices.reserve(nnz() + other.nnz());
  new_data.reserve(nnz() + other.nnz());

  for (int col = 0; col < n; ++col) {
    new_indptr[col] = static_cast<int>(new_indices.size());
    int a_start = a_indptr[col], a_end = a_indptr[col + 1];
    int b_start = b_indptr[col], b_end = b_indptr[col + 1];
    int ai = a_start, bi = b_start;

    while (ai < a_end && bi < b_end) {
      if (a_indices[ai] < b_indices[bi]) {
        new_indices.push_back(a_indices[ai]);
        new_data.push_back(data_[ai]);
        ++ai;
      } else if (a_indices[ai] > b_indices[bi]) {
        new_indices.push_back(b_indices[bi]);
        new_data.push_back(other.data_[bi]);
        ++bi;
      } else {
        // Same row: sum values, drop if zero
        double sum = data_[ai] + other.data_[bi];
        if (sum != 0.0) {
          new_indices.push_back(a_indices[ai]);
          new_data.push_back(sum);
        }
        ++ai;
        ++bi;
      }
    }
    while (ai < a_end) {
      new_indices.push_back(a_indices[ai]);
      new_data.push_back(data_[ai]);
      ++ai;
    }
    while (bi < b_end) {
      new_indices.push_back(b_indices[bi]);
      new_data.push_back(other.data_[bi]);
      ++bi;
    }
  }
  new_indptr[n] = static_cast<int>(new_indices.size());

  CSCPattern pattern(n, std::move(new_indptr), std::move(new_indices));
  return CSCMatrix(std::move(pattern), std::move(new_data));
}

CSCMatrix& CSCMatrix::operator+=(const CSCMatrix& other)
{
  int n = pattern_.n();
  if (other.n() != n) {
    fatal_error(fmt::format(
      "Cannot add CSC matrices with different dimensions ({} vs {})", n,
      other.n()));
  }

  // Fast path: if patterns are identical, add values in-place.
  if (pattern_ == other.pattern()) {
    bool has_zero = false;
    for (int k = 0; k < nnz(); ++k) {
      data_[k] += other.data_[k];
      has_zero |= (data_[k] == 0.0);
    }
    if (!has_zero)
      return *this;

    // Cancellation produced zeros — rebuild to drop them
    const auto& ip = indptr();
    const auto& ix = indices();
    vector<int> new_indptr(n + 1);
    vector<int> new_indices;
    vector<double> new_data;
    new_indices.reserve(nnz());
    new_data.reserve(nnz());
    for (int col = 0; col < n; ++col) {
      new_indptr[col] = static_cast<int>(new_indices.size());
      for (int p = ip[col]; p < ip[col + 1]; ++p) {
        if (data_[p] != 0.0) {
          new_indices.push_back(ix[p]);
          new_data.push_back(data_[p]);
        }
      }
    }
    new_indptr[n] = static_cast<int>(new_indices.size());
    pattern_ = CSCPattern(n, std::move(new_indptr), std::move(new_indices));
    data_ = std::move(new_data);
    return *this;
  }

  // General path: fall back to operator+ (which drops zeros)
  *this = *this + other;
  return *this;
}

CSCMatrix operator*(double scalar, const CSCMatrix& mat)
{
  CSCMatrix result(mat);
  result.scale(scalar);
  return result;
}

void CSCMatrix::scale(double scalar)
{
  for (size_t k = 0; k < data_.size(); ++k) {
    data_[k] *= scalar;
  }
}

} // namespace openmc
