#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include "openmc/sparse_matrix.h"

using namespace openmc;

TEST_CASE("CSCPattern construction and accessors")
{
  // 3x3 pattern:
  //   col 0: rows 0, 2
  //   col 1: row 1
  //   col 2: rows 0, 2
  vector<int> indptr = {0, 2, 3, 5};
  vector<int> indices = {0, 2, 1, 0, 2};
  CSCPattern pat(3, indptr, indices);

  CHECK(pat.n() == 3);
  CHECK(pat.nnz() == 5);
  CHECK(pat.indptr() == indptr);
  CHECK(pat.indices() == indices);
}

TEST_CASE("CSCPattern default constructor")
{
  CSCPattern pat;
  CHECK(pat.n() == 0);
  CHECK(pat.nnz() == 0);
}

TEST_CASE("CSCPattern operator==")
{
  vector<int> indptr = {0, 2, 3, 5};
  vector<int> indices = {0, 2, 1, 0, 2};
  CSCPattern a(3, indptr, indices);
  CSCPattern b(3, indptr, indices);

  CHECK(a == b);
  CHECK_FALSE(a != b);

  // Different dimension
  CSCPattern c(4, {0, 0, 0, 0, 0}, {});
  CHECK(a != c);

  // Same dimension but different pattern
  CSCPattern d(3, {0, 1, 2, 3}, {0, 1, 2});
  CHECK(a != d);
}

TEST_CASE("CSCPattern with_diagonal inserts missing diagonals")
{
  // Pattern with no diagonal entries:
  //   col 0: row 2
  //   col 1: (empty)
  //   col 2: row 0
  CSCPattern pat(3, {0, 1, 1, 2}, {2, 0});
  auto wd = pat.with_diagonal();

  CHECK(wd.n() == 3);
  // Original 2 entries + 3 missing diagonals = 5
  CHECK(wd.nnz() == 5);

  // Verify diagonals are present and rows stay sorted
  const auto& ip = wd.indptr();
  const auto& ix = wd.indices();

  // col 0: rows 0(new diag), 2
  CHECK(ip[0] == 0);
  CHECK(ix[0] == 0);
  CHECK(ix[1] == 2);
  // col 1: row 1(new diag)
  CHECK(ip[1] == 2);
  CHECK(ix[2] == 1);
  // col 2: rows 0, 2(new diag)
  CHECK(ip[2] == 3);
  CHECK(ix[3] == 0);
  CHECK(ix[4] == 2);
  CHECK(ip[3] == 5);
}

TEST_CASE("CSCPattern with_diagonal preserves existing diagonals")
{
  // Full diagonal already present
  CSCPattern pat(3, {0, 2, 4, 6}, {0, 1, 1, 2, 0, 2});
  auto wd = pat.with_diagonal();

  CHECK(wd.n() == 3);
  CHECK(wd.nnz() == 6); // unchanged
  CHECK(wd == pat);
}

TEST_CASE("CSCPattern with_diagonal partial diagonals")
{
  // 3x3 pattern: col 0 has diag, col 1 missing, col 2 has diag
  //   col 0: rows 0, 2
  //   col 1: row 0
  //   col 2: rows 1, 2
  CSCPattern pat(3, {0, 2, 3, 5}, {0, 2, 0, 1, 2});
  auto wd = pat.with_diagonal();

  CHECK(wd.nnz() == 6); // 5 + 1 missing diagonal for col 1

  // col 1 should now have rows 0, 1 (diagonal inserted in sorted order)
  const auto& ip = wd.indptr();
  const auto& ix = wd.indices();
  CHECK(ix[ip[1]] == 0);
  CHECK(ix[ip[1] + 1] == 1);
}

TEST_CASE("CSCMatrix construction and accessors")
{
  // 3x3 matrix:
  //   [1  0  2]
  //   [0  3  0]
  //   [4  0  5]
  vector<int> indptr = {0, 2, 3, 5};
  vector<int> indices = {0, 2, 1, 0, 2};
  vector<double> data = {1.0, 4.0, 3.0, 2.0, 5.0};
  CSCMatrix mat(3, indptr, indices, data);

  CHECK(mat.n() == 3);
  CHECK(mat.nnz() == 5);
  CHECK(mat.indptr() == indptr);
  CHECK(mat.indices() == indices);
  CHECK(mat.data() == data);

  // Verify pattern matches equivalent standalone CSCPattern
  CSCPattern pat(3, indptr, indices);
  CHECK(mat.pattern() == pat);
}

TEST_CASE("CSCMatrix default constructor")
{
  CSCMatrix mat;
  CHECK(mat.n() == 0);
  CHECK(mat.nnz() == 0);
}

TEST_CASE("CSCPattern rejects malformed inputs")
{
  // indptr size mismatch (n=3 requires 4 entries)
  CHECK_THROWS_AS(CSCPattern(3, {0, 1, 2}, {0, 1}), std::invalid_argument);

  // indptr[0] must be 0
  CHECK_THROWS_AS(CSCPattern(2, {1, 1, 2}, {0, 1}), std::invalid_argument);

  // indptr[n] must equal nnz
  CHECK_THROWS_AS(CSCPattern(2, {0, 1, 3}, {0, 1}), std::invalid_argument);

  // Row index out of bounds
  CHECK_THROWS_AS(CSCPattern(2, {0, 1, 2}, {0, 5}), std::invalid_argument);

  // Unsorted row indices within a column
  CHECK_THROWS_AS(CSCPattern(3, {0, 2, 2, 2}, {2, 0}), std::invalid_argument);

  // Duplicate row indices within a column (not strictly ascending)
  CHECK_THROWS_AS(CSCPattern(3, {0, 2, 2, 2}, {1, 1}), std::invalid_argument);
}

TEST_CASE("CSCMatrix rejects data/nnz size mismatch")
{
  // Pattern has 2 nonzeros, only 1 data value supplied
  CHECK_THROWS_AS(
    CSCMatrix(2, {0, 1, 2}, {0, 1}, {1.0}), std::invalid_argument);
}
