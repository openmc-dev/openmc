#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/tensor.h"

using namespace openmc;
using namespace openmc::tensor;

// ============================================================================
// Tensor constructors
// ============================================================================

TEST_CASE("Tensor default constructor")
{
  Tensor<double> t;
  REQUIRE(t.size() == 0);
  REQUIRE(t.empty());
  REQUIRE(t.shape().empty());
}

TEST_CASE("Tensor shape constructor")
{
  Tensor<double> t1({5});
  REQUIRE(t1.size() == 5);
  REQUIRE(t1.shape().size() == 1);
  REQUIRE(t1.shape(0) == 5);

  Tensor<double> t2({3, 4});
  REQUIRE(t2.size() == 12);
  REQUIRE(t2.shape().size() == 2);
  REQUIRE(t2.shape(0) == 3);
  REQUIRE(t2.shape(1) == 4);

  Tensor<int> t3({2, 3, 4});
  REQUIRE(t3.size() == 24);
  REQUIRE(t3.shape().size() == 3);
}

TEST_CASE("Tensor shape + fill constructor")
{
  Tensor<double> t({2, 3}, 7.0);
  REQUIRE(t.size() == 6);
  for (size_t i = 0; i < t.size(); ++i)
    REQUIRE(t[i] == 7.0);
}

TEST_CASE("Tensor pointer constructor")
{
  double vals[] = {1.0, 2.0, 3.0, 4.0};
  Tensor<double> t(vals, 4);
  REQUIRE(t.size() == 4);
  REQUIRE(t.shape(0) == 4);
  REQUIRE(t[0] == 1.0);
  REQUIRE(t[1] == 2.0);
  REQUIRE(t[2] == 3.0);
  REQUIRE(t[3] == 4.0);
}

TEST_CASE("Tensor copy and move")
{
  Tensor<double> a({2, 3}, 5.0);
  Tensor<double> b(a);
  REQUIRE(b.size() == 6);
  REQUIRE(b(0, 0) == 5.0);
  // Modifying copy doesn't affect original
  b(0, 0) = 99.0;
  REQUIRE(a(0, 0) == 5.0);

  Tensor<double> c(std::move(b));
  REQUIRE(c(0, 0) == 99.0);
  REQUIRE(c.size() == 6);
}

// ============================================================================
// Tensor indexing
// ============================================================================

TEST_CASE("Tensor 1D indexing")
{
  Tensor<int> t({4}, 0);
  t[0] = 10;
  t[1] = 20;
  t[2] = 30;
  t[3] = 40;
  REQUIRE(t(0) == 10);
  REQUIRE(t(1) == 20);
  REQUIRE(t(2) == 30);
  REQUIRE(t(3) == 40);
}

TEST_CASE("Tensor 2D indexing (row-major)")
{
  // Layout: [[1, 2, 3], [4, 5, 6]]
  Tensor<int> t({2, 3}, 0);
  int val = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      t(i, j) = val++;

  REQUIRE(t(0, 0) == 1);
  REQUIRE(t(0, 2) == 3);
  REQUIRE(t(1, 0) == 4);
  REQUIRE(t(1, 2) == 6);
  // Flat index should match row-major order
  REQUIRE(t[0] == 1);
  REQUIRE(t[3] == 4);
  REQUIRE(t[5] == 6);
}

TEST_CASE("Tensor 3D indexing")
{
  // 2x3x4 tensor
  Tensor<int> t({2, 3, 4}, 0);
  t(1, 2, 3) = 42;
  // Flat index: 1*12 + 2*4 + 3 = 23
  REQUIRE(t[23] == 42);
  REQUIRE(t(1, 2, 3) == 42);
}

// ============================================================================
// Tensor assignment
// ============================================================================

TEST_CASE("Tensor initializer_list assignment")
{
  Tensor<double> t;
  t = {1.0, 2.0, 3.0};
  REQUIRE(t.size() == 3);
  REQUIRE(t.shape(0) == 3);
  REQUIRE(t[0] == 1.0);
  REQUIRE(t[2] == 3.0);
}

// ============================================================================
// Tensor mutation
// ============================================================================

TEST_CASE("Tensor resize")
{
  Tensor<double> t({2, 3}, 1.0);
  REQUIRE(t.size() == 6);
  t.resize({4, 5});
  REQUIRE(t.size() == 20);
  REQUIRE(t.shape(0) == 4);
  REQUIRE(t.shape(1) == 5);
}

TEST_CASE("Tensor reshape")
{
  Tensor<int> t({12}, 0);
  for (size_t i = 0; i < 12; ++i)
    t[i] = static_cast<int>(i);

  t.reshape({3, 4});
  REQUIRE(t.shape(0) == 3);
  REQUIRE(t.shape(1) == 4);
  REQUIRE(t.size() == 12);
  // Data unchanged, just reinterpreted
  REQUIRE(t(0, 0) == 0);
  REQUIRE(t(1, 0) == 4);  // row 1, col 0 = flat index 4
  REQUIRE(t(2, 3) == 11); // row 2, col 3 = flat index 11
}

TEST_CASE("Tensor fill")
{
  Tensor<double> t({3, 3}, 0.0);
  t.fill(42.0);
  for (size_t i = 0; i < t.size(); ++i)
    REQUIRE(t[i] == 42.0);
}

// ============================================================================
// Tensor iterators
// ============================================================================

TEST_CASE("Tensor iterators")
{
  Tensor<int> t({4}, 0);
  t = {10, 20, 30, 40};
  int sum = 0;
  for (auto val : t)
    sum += val;
  REQUIRE(sum == 100);
}

// ============================================================================
// Tensor reductions
// ============================================================================

TEST_CASE("Tensor sum (full)")
{
  Tensor<double> t({3}, 0.0);
  t = {1.0, 2.0, 3.0};
  REQUIRE(t.sum() == 6.0);
}

TEST_CASE("Tensor sum (axis) on 2D")
{
  // [[1, 2, 3],
  //  [4, 5, 6]]
  Tensor<int> t({2, 3}, 0);
  int v = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      t(i, j) = v++;

  // Sum along axis 0 -> [5, 7, 9]
  Tensor<int> s0 = t.sum(0);
  REQUIRE(s0.size() == 3);
  REQUIRE(s0[0] == 5);
  REQUIRE(s0[1] == 7);
  REQUIRE(s0[2] == 9);

  // Sum along axis 1 -> [6, 15]
  Tensor<int> s1 = t.sum(1);
  REQUIRE(s1.size() == 2);
  REQUIRE(s1[0] == 6);
  REQUIRE(s1[1] == 15);
}

TEST_CASE("Tensor sum (axis) on 3D")
{
  // 2x3x2 tensor filled with sequential values 1..12
  Tensor<int> t({2, 3, 2}, 0);
  int v = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 2; ++k)
        t(i, j, k) = v++;

  // Sum along axis 1 (middle) -> 2x2, each sums 3 values
  // [0,0]: t(0,0,0)+t(0,1,0)+t(0,2,0) = 1+3+5 = 9
  // [0,1]: t(0,0,1)+t(0,1,1)+t(0,2,1) = 2+4+6 = 12
  // [1,0]: t(1,0,0)+t(1,1,0)+t(1,2,0) = 7+9+11 = 27
  // [1,1]: t(1,0,1)+t(1,1,1)+t(1,2,1) = 8+10+12 = 30
  Tensor<int> s = t.sum(1);
  REQUIRE(s.shape(0) == 2);
  REQUIRE(s.shape(1) == 2);
  REQUIRE(s(0, 0) == 9);
  REQUIRE(s(0, 1) == 12);
  REQUIRE(s(1, 0) == 27);
  REQUIRE(s(1, 1) == 30);
}

TEST_CASE("Tensor prod")
{
  Tensor<int> t({4}, 0);
  t = {1, 2, 3, 4};
  REQUIRE(t.prod() == 24);
}

TEST_CASE("Tensor any and all")
{
  Tensor<bool> t({4}, false);
  REQUIRE(!t.any());
  REQUIRE(!t.all());

  // Set one element true
  t.data()[0] = true;
  REQUIRE(t.any());
  REQUIRE(!t.all());

  // Set all true
  for (size_t i = 0; i < t.size(); ++i)
    t.data()[i] = true;
  REQUIRE(t.any());
  REQUIRE(t.all());
}

TEST_CASE("Tensor argmin")
{
  Tensor<double> t({5}, 0.0);
  t = {3.0, 1.0, 4.0, 0.5, 2.0};
  REQUIRE(t.argmin() == 3);
}

TEST_CASE("Tensor flip")
{
  Tensor<int> t({5}, 0);
  t = {1, 2, 3, 4, 5};
  Tensor<int> f = t.flip(0);
  REQUIRE(f[0] == 5);
  REQUIRE(f[1] == 4);
  REQUIRE(f[2] == 3);
  REQUIRE(f[3] == 2);
  REQUIRE(f[4] == 1);
}

TEST_CASE("Tensor flip 2D")
{
  // [[1, 2], [3, 4], [5, 6]]
  Tensor<int> t({3, 2}, 0);
  t(0, 0) = 1;
  t(0, 1) = 2;
  t(1, 0) = 3;
  t(1, 1) = 4;
  t(2, 0) = 5;
  t(2, 1) = 6;

  // Flip axis 0 reverses rows -> [[5,6],[3,4],[1,2]]
  Tensor<int> f = t.flip(0);
  REQUIRE(f(0, 0) == 5);
  REQUIRE(f(0, 1) == 6);
  REQUIRE(f(1, 0) == 3);
  REQUIRE(f(2, 0) == 1);
}

// ============================================================================
// Tensor operators
// ============================================================================

TEST_CASE("Tensor scalar compound assignment")
{
  Tensor<double> t({3}, 0.0);
  t = {2.0, 4.0, 6.0};

  t += 1.0;
  REQUIRE(t[0] == 3.0);
  REQUIRE(t[1] == 5.0);

  t -= 1.0;
  REQUIRE(t[0] == 2.0);

  t *= 3.0;
  REQUIRE(t[0] == 6.0);
  REQUIRE(t[1] == 12.0);

  t /= 2.0;
  REQUIRE(t[0] == 3.0);
  REQUIRE(t[1] == 6.0);
}

TEST_CASE("Tensor element-wise arithmetic")
{
  Tensor<double> a({3}, 0.0);
  Tensor<double> b({3}, 0.0);
  a = {1.0, 2.0, 3.0};
  b = {4.0, 5.0, 6.0};

  Tensor<double> c = a + b;
  REQUIRE(c[0] == 5.0);
  REQUIRE(c[1] == 7.0);
  REQUIRE(c[2] == 9.0);

  c = a - b;
  REQUIRE(c[0] == -3.0);

  c = a / b;
  REQUIRE(c[0] == 0.25);

  c = a * b;
  REQUIRE(c[0] == 4.0);
}

TEST_CASE("Tensor scalar arithmetic")
{
  Tensor<double> a({3}, 0.0);
  a = {1.0, 2.0, 3.0};

  Tensor<double> b = a + 10.0;
  REQUIRE(b[0] == 11.0);
  REQUIRE(b[2] == 13.0);

  b = a - 1.0;
  REQUIRE(b[0] == 0.0);

  b = a * 2.0;
  REQUIRE(b[0] == 2.0);
  REQUIRE(b[2] == 6.0);

  // Non-member scalar * tensor (commutativity)
  b = 2.0 * a;
  REQUIRE(b[0] == 2.0);
  REQUIRE(b[2] == 6.0);

  // Non-member scalar + tensor
  b = 10.0 + a;
  REQUIRE(b[0] == 11.0);
}

TEST_CASE("Tensor compound arithmetic with tensor")
{
  Tensor<double> a({3}, 0.0);
  Tensor<double> b({3}, 0.0);
  a = {1.0, 2.0, 3.0};
  b = {10.0, 20.0, 30.0};
  a += b;
  REQUIRE(a[0] == 11.0);
  REQUIRE(a[1] == 22.0);
  REQUIRE(a[2] == 33.0);

  a = {1.0, 2.0, 3.0};
  b -= a;
  REQUIRE(b[0] == 9.0);
  REQUIRE(b[1] == 18.0);
  REQUIRE(b[2] == 27.0);

  b = {10.0, 20.0, 30.0};
  a *= b;
  REQUIRE(a[0] == 10.0);
  REQUIRE(a[1] == 40.0);
  REQUIRE(a[2] == 90.0);

  a = {1.0, 2.0, 3.0};
  b /= a;
  REQUIRE(b[0] == 10.0);
  REQUIRE(b[1] == 10.0);
  REQUIRE(b[2] == 10.0);
}

TEST_CASE("Tensor comparison operators")
{
  Tensor<double> t({4}, 0.0);
  t = {1.0, 2.0, 3.0, 4.0};

  Tensor<bool> r = t < 3.0;
  REQUIRE(r.data()[0] == true);
  REQUIRE(r.data()[1] == true);
  REQUIRE(r.data()[2] == false);
  REQUIRE(r.data()[3] == false);

  r = t >= 3.0;
  REQUIRE(r.data()[0] == false);
  REQUIRE(r.data()[2] == true);
  REQUIRE(r.data()[3] == true);

  r = t <= 2.0;
  REQUIRE(r.data()[0] == true);
  REQUIRE(r.data()[1] == true);
  REQUIRE(r.data()[2] == false);

  r = t > 3.0;
  REQUIRE(r.data()[0] == false);
  REQUIRE(r.data()[3] == true);
}

TEST_CASE("Tensor element-wise comparison")
{
  Tensor<double> a({3}, 0.0);
  Tensor<double> b({3}, 0.0);
  a = {1.0, 5.0, 3.0};
  b = {2.0, 4.0, 3.0};

  Tensor<bool> r = a < b;
  REQUIRE(r.data()[0] == true);
  REQUIRE(r.data()[1] == false);
  REQUIRE(r.data()[2] == false);
}

TEST_CASE("Tensor mixed-type multiply")
{
  Tensor<int> a({3}, 0);
  Tensor<double> b({3}, 0.0);
  a = {2, 3, 4};
  b = {1.5, 2.5, 3.5};

  Tensor<double> c = a * b;
  REQUIRE(c[0] == 3.0);
  REQUIRE(c[1] == 7.5);
  REQUIRE(c[2] == 14.0);
}

TEST_CASE("Tensor mixed-type divide")
{
  Tensor<double> a({3}, 0.0);
  Tensor<int> b({3}, 0);
  a = {10.0, 20.0, 30.0};
  b = {2, 4, 5};

  Tensor<double> c = a / b;
  REQUIRE(c[0] == 5.0);
  REQUIRE(c[1] == 5.0);
  REQUIRE(c[2] == 6.0);
}

// ============================================================================
// Tensor bool specialization
// ============================================================================

TEST_CASE("Tensor<bool> storage")
{
  // Tensor<bool> uses unsigned char internally to avoid std::vector<bool> proxy
  Tensor<bool> t({4}, false);
  t.data()[0] = true;
  t.data()[2] = true;
  REQUIRE(t.any());
  REQUIRE(!t.all());
  REQUIRE(t.data()[0] == true);
  REQUIRE(t.data()[1] == false);
}

// ============================================================================
// View (via Tensor accessors)
// ============================================================================

TEST_CASE("Tensor slice axis 0 (2D)")
{
  // [[1, 2, 3], [4, 5, 6]]
  Tensor<int> t({2, 3}, 0);
  int v = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      t(i, j) = v++;

  auto r0 = t.slice(0);
  REQUIRE(r0.size() == 3);
  REQUIRE(r0[0] == 1);
  REQUIRE(r0[1] == 2);
  REQUIRE(r0[2] == 3);

  auto r1 = t.slice(1);
  REQUIRE(r1[0] == 4);
  REQUIRE(r1[1] == 5);
  REQUIRE(r1[2] == 6);

  // Writing through view modifies the tensor
  r0[1] = 99;
  REQUIRE(t(0, 1) == 99);
}

TEST_CASE("Tensor slice axis 1 (2D)")
{
  // [[1, 2], [3, 4], [5, 6]]
  Tensor<int> t({3, 2}, 0);
  t(0, 0) = 1;
  t(0, 1) = 2;
  t(1, 0) = 3;
  t(1, 1) = 4;
  t(2, 0) = 5;
  t(2, 1) = 6;

  auto c0 = t.slice(all, 0);
  REQUIRE(c0.size() == 3);
  REQUIRE(c0[0] == 1);
  REQUIRE(c0[1] == 3);
  REQUIRE(c0[2] == 5);

  auto c1 = t.slice(all, 1);
  REQUIRE(c1[0] == 2);
  REQUIRE(c1[1] == 4);
  REQUIRE(c1[2] == 6);

  // Write through column view
  c1[0] = 77;
  REQUIRE(t(0, 1) == 77);
}

TEST_CASE("Tensor slice with range")
{
  Tensor<int> t({6}, 0);
  t = {10, 20, 30, 40, 50, 60};

  // range(start, end)
  auto s = t.slice(range(1, 4));
  REQUIRE(s.size() == 3);
  REQUIRE(s[0] == 20);
  REQUIRE(s[1] == 30);
  REQUIRE(s[2] == 40);

  // range(end) from start — range(3) means [0, 3)
  auto s2 = t.slice(range(3));
  REQUIRE(s2.size() == 3);
  REQUIRE(s2[0] == 10);
  REQUIRE(s2[2] == 30);

  // range(start, SIZE_MAX) to end
  auto s3 = t.slice(range(3, 6));
  REQUIRE(s3.size() == 3);
  REQUIRE(s3[0] == 40);
  REQUIRE(s3[2] == 60);

  // Write through slice
  s[0] = 99;
  REQUIRE(t[1] == 99);
}

TEST_CASE("Tensor flat view")
{
  Tensor<int> t({2, 3}, 0);
  int v = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      t(i, j) = v++;

  auto f = t.flat();
  REQUIRE(f.size() == 6);
  REQUIRE(f[0] == 1);
  REQUIRE(f[5] == 6);
}

TEST_CASE("Tensor slice on 3D")
{
  // 2x3x4 tensor
  Tensor<int> t({2, 3, 4}, 0);
  int v = 0;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 4; ++k)
        t(i, j, k) = v++;

  // slice(1) -> fix axis 0 at 1 -> 3x4 view
  auto s = t.slice(1);
  REQUIRE(s.size() == 12);
  // t(1,0,0) = 12, t(1,0,1) = 13, ...
  REQUIRE(s(0, 0) == 12);
  REQUIRE(s(0, 1) == 13);
  REQUIRE(s(2, 3) == 23);

  // slice(all, 2) -> fix axis 1 at 2 -> 2x4 view
  auto s2 = t.slice(all, 2);
  REQUIRE(s2.size() == 8);
  // t(0,2,0)=8, t(0,2,1)=9, t(1,2,0)=20
  REQUIRE(s2(0, 0) == 8);
  REQUIRE(s2(0, 1) == 9);
  REQUIRE(s2(1, 0) == 20);
}

TEST_CASE("Tensor multi-axis slice")
{
  // 2x3x4 tensor with sequential values
  Tensor<int> t({2, 3, 4}, 0);
  int v = 0;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 4; ++k)
        t(i, j, k) = v++;

  // slice(1, 2) -> fix axes 0 and 1 -> 1D view of 4 elements
  // Equivalent to numpy t[1, 2, :] -> t(1,2,0..3) = [20, 21, 22, 23]
  auto s = t.slice(1, 2);
  REQUIRE(s.size() == 4);
  REQUIRE(s[0] == 20);
  REQUIRE(s[1] == 21);
  REQUIRE(s[3] == 23);

  // slice(all, 1, range(1, 3)) -> keep axis 0, fix axis 1, range on axis 2
  // Equivalent to numpy t[:, 1, 1:3]
  // t(0,1,1)=5, t(0,1,2)=6, t(1,1,1)=17, t(1,1,2)=18
  auto s2 = t.slice(all, 1, range(1, 3));
  REQUIRE(s2.size() == 4);
  REQUIRE(s2.ndim() == 2);
  REQUIRE(s2(0, 0) == 5);
  REQUIRE(s2(0, 1) == 6);
  REQUIRE(s2(1, 0) == 17);
  REQUIRE(s2(1, 1) == 18);

  // slice(0, range(0, 2)) -> fix axis 0 at 0, range on axis 1
  // Equivalent to numpy t[0, 0:2, :] -> shape (2, 4)
  auto s3 = t.slice(0, range(0, 2));
  REQUIRE(s3.ndim() == 2);
  REQUIRE(s3.shape(0) == 2);
  REQUIRE(s3.shape(1) == 4);
  REQUIRE(s3(0, 0) == 0); // t(0,0,0)
  REQUIRE(s3(1, 3) == 7); // t(0,1,3)
}

// ============================================================================
// View assignment and arithmetic
// ============================================================================

TEST_CASE("View scalar assignment (fill)")
{
  Tensor<double> t({2, 3}, 0.0);
  auto r = t.slice(0);
  r = 7.0;
  REQUIRE(t(0, 0) == 7.0);
  REQUIRE(t(0, 1) == 7.0);
  REQUIRE(t(0, 2) == 7.0);
  REQUIRE(t(1, 0) == 0.0); // Other row unchanged
}

TEST_CASE("View initializer_list assignment")
{
  Tensor<double> t({2, 3}, 0.0);
  auto r = t.slice(1);
  r = {10.0, 20.0, 30.0};
  REQUIRE(t(1, 0) == 10.0);
  REQUIRE(t(1, 1) == 20.0);
  REQUIRE(t(1, 2) == 30.0);
}

TEST_CASE("View copy assignment (deep copy)")
{
  Tensor<double> t({2, 3}, 0.0);
  t.slice(0) = {1.0, 2.0, 3.0};
  t.slice(1) = {4.0, 5.0, 6.0};

  // Copy row 0 into row 1
  t.slice(1) = t.slice(0);
  REQUIRE(t(1, 0) == 1.0);
  REQUIRE(t(1, 1) == 2.0);
  REQUIRE(t(1, 2) == 3.0);
}

TEST_CASE("View compound operators")
{
  Tensor<double> t({2, 3}, 0.0);
  t.slice(0) = {1.0, 2.0, 3.0};

  t.slice(0) *= 2.0;
  REQUIRE(t(0, 0) == 2.0);
  REQUIRE(t(0, 1) == 4.0);

  t.slice(0) /= 2.0;
  REQUIRE(t(0, 0) == 1.0);
  REQUIRE(t(0, 1) == 2.0);
}

TEST_CASE("View assignment from tensor")
{
  Tensor<double> t({2, 3}, 0.0);
  Tensor<double> vals({3}, 0.0);
  vals = {7.0, 8.0, 9.0};

  t.slice(1) = vals;
  REQUIRE(t(1, 0) == 7.0);
  REQUIRE(t(1, 1) == 8.0);
  REQUIRE(t(1, 2) == 9.0);
}

TEST_CASE("View compound addition from tensor")
{
  Tensor<double> t({2, 3}, 0.0);
  t.slice(0) = {1.0, 2.0, 3.0};
  Tensor<double> vals({3}, 0.0);
  vals = {10.0, 20.0, 30.0};

  t.slice(0) += vals;
  REQUIRE(t(0, 0) == 11.0);
  REQUIRE(t(0, 1) == 22.0);
  REQUIRE(t(0, 2) == 33.0);
}

TEST_CASE("View sum")
{
  Tensor<double> t({2, 3}, 0.0);
  t.slice(0) = {1.0, 2.0, 3.0};
  t.slice(1) = {4.0, 5.0, 6.0};

  REQUIRE(t.slice(0).sum() == 6.0);
  REQUIRE(t.slice(1).sum() == 15.0);
}

TEST_CASE("View iteration")
{
  Tensor<int> t({2, 3}, 0);
  t.slice(0) = {1, 2, 3};

  int sum = 0;
  for (auto val : t.slice(0))
    sum += val;
  REQUIRE(sum == 6);
}

TEST_CASE("View sub-slice")
{
  Tensor<int> t({6}, 0);
  t = {10, 20, 30, 40, 50, 60};

  auto s = t.slice(range(1, 5));  // [20, 30, 40, 50]
  auto ss = s.slice(range(1, 3)); // [30, 40]
  REQUIRE(ss.size() == 2);
  REQUIRE(ss[0] == 30);
  REQUIRE(ss[1] == 40);
}

TEST_CASE("Tensor from View")
{
  Tensor<double> t({2, 3}, 0.0);
  t.slice(0) = {1.0, 2.0, 3.0};

  // Construct a new tensor from a view (copies data)
  Tensor<double> t2(t.slice(0));
  REQUIRE(t2.size() == 3);
  REQUIRE(t2[0] == 1.0);
  REQUIRE(t2[2] == 3.0);

  // Modifying the new tensor doesn't affect the original
  t2[0] = 99.0;
  REQUIRE(t(0, 0) == 1.0);
}

// ============================================================================
// Const View
// ============================================================================

TEST_CASE("Const tensor produces const views")
{
  Tensor<double> t({2, 3}, 0.0);
  int v = 1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 3; ++j)
      t(i, j) = v++;

  const Tensor<double>& ct = t;
  auto r = ct.slice(0); // View<const double>
  REQUIRE(r[0] == 1.0);
  REQUIRE(r[2] == 3.0);

  auto c = ct.slice(all, 1);
  REQUIRE(c[0] == 2.0);
  REQUIRE(c[1] == 5.0);
}

// ============================================================================
// StaticTensor2D
// ============================================================================

TEST_CASE("StaticTensor2D basics")
{
  StaticTensor2D<double, 3, 4> t;
  REQUIRE(t.size() == 12);
  REQUIRE(t.shape()[0] == 3);
  REQUIRE(t.shape()[1] == 4);

  // Default-initialized to zero
  REQUIRE(t(0, 0) == 0.0);

  t(1, 2) = 42.0;
  REQUIRE(t(1, 2) == 42.0);
  // Flat data: row 1, col 2 = index 1*4 + 2 = 6
  REQUIRE(t.data()[6] == 42.0);
}

TEST_CASE("StaticTensor2D fill")
{
  StaticTensor2D<int, 2, 3> t;
  t.fill(5);
  for (size_t i = 0; i < t.size(); ++i)
    REQUIRE(t.data()[i] == 5);
}

TEST_CASE("StaticTensor2D iteration")
{
  StaticTensor2D<int, 2, 3> t;
  t.fill(1);
  int sum = 0;
  for (auto val : t)
    sum += val;
  REQUIRE(sum == 6);
}

TEST_CASE("StaticTensor2D slice")
{
  StaticTensor2D<int, 3, 2> t;
  t(0, 0) = 1;
  t(0, 1) = 2;
  t(1, 0) = 3;
  t(1, 1) = 4;
  t(2, 0) = 5;
  t(2, 1) = 6;

  // slice(1) = row 1 (fix axis 0 at 1)
  auto r1 = t.slice(1);
  REQUIRE(r1.size() == 2);
  REQUIRE(r1[0] == 3);
  REQUIRE(r1[1] == 4);

  // slice(all, 0) = column 0 (fix axis 1 at 0)
  auto c0 = t.slice(all, 0);
  REQUIRE(c0.size() == 3);
  REQUIRE(c0[0] == 1);
  REQUIRE(c0[1] == 3);
  REQUIRE(c0[2] == 5);
}

TEST_CASE("StaticTensor2D flat view")
{
  StaticTensor2D<double, 2, 2> t;
  t(0, 0) = 1.0;
  t(0, 1) = 2.0;
  t(1, 0) = 3.0;
  t(1, 1) = 4.0;

  auto f = t.flat();
  REQUIRE(f.size() == 4);
  f = 0.0;
  REQUIRE(t(0, 0) == 0.0);
  REQUIRE(t(1, 1) == 0.0);
}

// ============================================================================
// Non-member functions
// ============================================================================

TEST_CASE("zeros")
{
  auto t = zeros<double>({3, 4});
  REQUIRE(t.size() == 12);
  for (size_t i = 0; i < t.size(); ++i)
    REQUIRE(t[i] == 0.0);
}

TEST_CASE("zeros_like")
{
  Tensor<double> a({2, 5}, 7.0);
  auto b = zeros_like(a);
  REQUIRE(b.size() == 10);
  REQUIRE(b.shape(0) == 2);
  REQUIRE(b.shape(1) == 5);
  for (size_t i = 0; i < b.size(); ++i)
    REQUIRE(b[i] == 0.0);
}

TEST_CASE("full_like")
{
  Tensor<int> a({4}, 0);
  auto b = full_like(a, 42);
  REQUIRE(b.size() == 4);
  for (size_t i = 0; i < b.size(); ++i)
    REQUIRE(b[i] == 42);
}

TEST_CASE("linspace")
{
  auto t = linspace(0.0, 1.0, 5);
  REQUIRE(t.size() == 5);
  REQUIRE(t[0] == 0.0);
  REQUIRE(t[4] == 1.0);
  REQUIRE_THAT(t[1], Catch::Matchers::WithinRel(0.25, 1e-12));
  REQUIRE_THAT(t[2], Catch::Matchers::WithinRel(0.5, 1e-12));
  REQUIRE_THAT(t[3], Catch::Matchers::WithinRel(0.75, 1e-12));
}

TEST_CASE("concatenate")
{
  Tensor<int> a({3}, 0);
  Tensor<int> b({2}, 0);
  a = {1, 2, 3};
  b = {4, 5};

  auto c = concatenate(a, b);
  REQUIRE(c.size() == 5);
  REQUIRE(c[0] == 1);
  REQUIRE(c[2] == 3);
  REQUIRE(c[3] == 4);
  REQUIRE(c[4] == 5);
}

TEST_CASE("log")
{
  Tensor<double> t({3}, 0.0);
  t = {1.0, std::exp(1.0), std::exp(2.0)};

  auto r = log(t);
  REQUIRE_THAT(r[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(r[1], Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(r[2], Catch::Matchers::WithinAbs(2.0, 1e-12));
}

TEST_CASE("abs")
{
  Tensor<double> t({4}, 0.0);
  t = {-3.0, -1.0, 0.0, 2.0};

  auto r = abs(t);
  REQUIRE(r[0] == 3.0);
  REQUIRE(r[1] == 1.0);
  REQUIRE(r[2] == 0.0);
  REQUIRE(r[3] == 2.0);
}

TEST_CASE("where")
{
  Tensor<bool> cond({4}, false);
  cond.data()[0] = true;
  cond.data()[2] = true;

  Tensor<double> vals({4}, 0.0);
  vals = {10.0, 20.0, 30.0, 40.0};

  auto r = where(cond, vals, -1.0);
  REQUIRE(r[0] == 10.0);
  REQUIRE(r[1] == -1.0);
  REQUIRE(r[2] == 30.0);
  REQUIRE(r[3] == -1.0);
}

TEST_CASE("nan_to_num")
{
  Tensor<double> t({4}, 0.0);
  t[0] = 1.0;
  t[1] = std::nan("");
  t[2] = std::numeric_limits<double>::infinity();
  t[3] = -std::numeric_limits<double>::infinity();

  auto r = nan_to_num(t);
  REQUIRE(r[0] == 1.0);
  REQUIRE(r[1] == 0.0);                                   // NaN -> 0
  REQUIRE(r[2] == std::numeric_limits<double>::max());    // +inf -> max
  REQUIRE(r[3] == std::numeric_limits<double>::lowest()); // -inf -> lowest
}

// ============================================================================
// is_tensor trait
// ============================================================================

TEST_CASE("is_tensor trait")
{
  REQUIRE(is_tensor<Tensor<double>>::value);
  REQUIRE(is_tensor<Tensor<int>>::value);
  REQUIRE(is_tensor<StaticTensor2D<double, 3, 3>>::value);
  REQUIRE(!is_tensor<double>::value);
  REQUIRE(!is_tensor<std::vector<double>>::value);
}
