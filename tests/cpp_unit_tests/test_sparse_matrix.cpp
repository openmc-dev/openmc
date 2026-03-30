#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/sparse_matrix.h"

using namespace openmc;
using Catch::Matchers::WithinAbs;

// Helper: build a simple 3x3 matrix
//   [1  0  2]
//   [0  3  0]
//   [4  0  5]
// In CSC (column-major):
//   col 0: rows 0,2 vals 1,4
//   col 1: row 1    val 3
//   col 2: rows 0,2 vals 2,5
static CSCMatrix make_A()
{
  return CSCMatrix::from_triplets(3, {0, 1, 2, 0, 2}, {0, 1, 0, 2, 2},
    {1.0, 3.0, 4.0, 2.0, 5.0});
}

// Helper: build a 3x3 matrix with different pattern
//   [0  6  0]
//   [7  0  0]
//   [0  0  8]
static CSCMatrix make_B()
{
  return CSCMatrix::from_triplets(3, {1, 0, 2}, {0, 1, 2}, {7.0, 6.0, 8.0});
}

TEST_CASE("CSCMatrix scalar multiply")
{
  auto A = make_A();

  auto C = 2.0 * A;
  REQUIRE(C.n() == 3);
  REQUIRE(C.nnz() == 5);

  // Check values: should be doubled
  const auto& d = C.data();
  CHECK_THAT(d[0], WithinAbs(2.0, 1e-15));  // A(0,0) * 2
  CHECK_THAT(d[1], WithinAbs(8.0, 1e-15));  // A(2,0) * 2
  CHECK_THAT(d[2], WithinAbs(6.0, 1e-15));  // A(1,1) * 2
  CHECK_THAT(d[3], WithinAbs(4.0, 1e-15));  // A(0,2) * 2
  CHECK_THAT(d[4], WithinAbs(10.0, 1e-15)); // A(2,2) * 2

  // Pattern must be identical
  CHECK(C.pattern() == A.pattern());
}

TEST_CASE("CSCMatrix scalar multiply by zero")
{
  auto A = make_A();
  auto C = 0.0 * A;
  REQUIRE(C.nnz() == 5); // structural zeros remain
  for (int k = 0; k < C.nnz(); ++k) {
    CHECK_THAT(C.data()[k], WithinAbs(0.0, 1e-15));
  }
}

TEST_CASE("CSCMatrix scale in-place")
{
  auto A = make_A();
  A.scale(3.0);
  const auto& d = A.data();
  CHECK_THAT(d[0], WithinAbs(3.0, 1e-15));
  CHECK_THAT(d[1], WithinAbs(12.0, 1e-15));
  CHECK_THAT(d[2], WithinAbs(9.0, 1e-15));
  CHECK_THAT(d[3], WithinAbs(6.0, 1e-15));
  CHECK_THAT(d[4], WithinAbs(15.0, 1e-15));
}

TEST_CASE("CSCMatrix operator+= same pattern")
{
  auto A = make_A();
  auto A2 = make_A();
  A += A2;

  // Same pattern: fast path, values doubled
  CHECK(A.nnz() == 5);
  const auto& d = A.data();
  CHECK_THAT(d[0], WithinAbs(2.0, 1e-15));
  CHECK_THAT(d[1], WithinAbs(8.0, 1e-15));
  CHECK_THAT(d[2], WithinAbs(6.0, 1e-15));
  CHECK_THAT(d[3], WithinAbs(4.0, 1e-15));
  CHECK_THAT(d[4], WithinAbs(10.0, 1e-15));
}

TEST_CASE("CSCMatrix operator+= different pattern")
{
  auto A = make_A();
  auto B = make_B();
  A += B;

  // Union of patterns: nnz should be 5 + 3 - 1 (overlap at (2,2)) = 7
  CHECK(A.n() == 3);
  CHECK(A.nnz() == 7);

  // Verify combined values by checking known entries
  // Reconstruct into dense for checking
  const auto& indptr = A.indptr();
  const auto& indices = A.indices();
  const auto& data = A.data();

  double dense[3][3] = {};
  for (int col = 0; col < 3; ++col) {
    for (int k = indptr[col]; k < indptr[col + 1]; ++k) {
      dense[indices[k]][col] = data[k];
    }
  }

  // A + B:
  // [1  6  2]
  // [7  3  0]
  // [4  0  13]
  CHECK_THAT(dense[0][0], WithinAbs(1.0, 1e-15));
  CHECK_THAT(dense[0][1], WithinAbs(6.0, 1e-15));
  CHECK_THAT(dense[0][2], WithinAbs(2.0, 1e-15));
  CHECK_THAT(dense[1][0], WithinAbs(7.0, 1e-15));
  CHECK_THAT(dense[1][1], WithinAbs(3.0, 1e-15));
  CHECK_THAT(dense[1][2], WithinAbs(0.0, 1e-15));
  CHECK_THAT(dense[2][0], WithinAbs(4.0, 1e-15));
  CHECK_THAT(dense[2][1], WithinAbs(0.0, 1e-15));
  CHECK_THAT(dense[2][2], WithinAbs(13.0, 1e-15));
}

TEST_CASE("CSCMatrix combined: A_decay + s * A_rxn pattern")
{
  // This is the actual use case: A = A_decay + s * A_rxn
  auto A_decay = make_A();
  auto A_rxn = make_B();
  double s = 1.5;

  auto A = A_decay + s * A_rxn;
  CHECK(A.n() == 3);

  // Reconstruct dense
  const auto& indptr = A.indptr();
  const auto& indices = A.indices();
  const auto& data = A.data();

  double dense[3][3] = {};
  for (int col = 0; col < 3; ++col) {
    for (int k = indptr[col]; k < indptr[col + 1]; ++k) {
      dense[indices[k]][col] = data[k];
    }
  }

  // A_decay + 1.5 * A_rxn:
  // [1+0    0+9    2+0  ]     [1    9    2  ]
  // [0+10.5 3+0    0+0  ]  =  [10.5 3    0  ]
  // [4+0    0+0    5+12 ]     [4    0    17  ]
  CHECK_THAT(dense[0][0], WithinAbs(1.0, 1e-15));
  CHECK_THAT(dense[0][1], WithinAbs(9.0, 1e-15));
  CHECK_THAT(dense[0][2], WithinAbs(2.0, 1e-15));
  CHECK_THAT(dense[1][0], WithinAbs(10.5, 1e-15));
  CHECK_THAT(dense[1][1], WithinAbs(3.0, 1e-15));
  CHECK_THAT(dense[2][0], WithinAbs(4.0, 1e-15));
  CHECK_THAT(dense[2][2], WithinAbs(17.0, 1e-15));
}

TEST_CASE("CSCMatrix empty matrix arithmetic")
{
  CSCMatrix empty;
  CHECK(empty.n() == 0);
  CHECK(empty.nnz() == 0);

  // scale empty
  empty.scale(5.0);
  CHECK(empty.nnz() == 0);

  // scalar * empty
  auto C = 3.0 * empty;
  CHECK(C.nnz() == 0);
}
