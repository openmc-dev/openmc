// Tests for the CORE-MATH backend of openmc/math.h.
//
// CORE-MATH's functions are correctly rounded, so for every input there is
// exactly one right answer and it must be produced on every platform and
// compiler. The known-answer vectors in portable_math_kat.h record the bits
// returned on the reference platform (Linux x86-64, GCC) for the special
// values and a sample of ordinary inputs; they are also the correctly rounded
// values, as checked against MPFR by tools/dev/portable_math_report.py. A
// platform that changes a single bit fails here before it can change a
// regression test result. The vectors are regenerated with
//
//   tools/dev/portable_math_report.py samples.txt --kat <this header>
//
// whenever the pinned CORE-MATH commit changes.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "openmc/math.h"
#include "portable_math_kat.h"

#ifndef OPENMC_COREMATH
#define OPENMC_PORTABLE_MATH_AVAILABLE 0
#else
#define OPENMC_PORTABLE_MATH_AVAILABLE 1
#endif

namespace {

double from_bits(uint64_t u)
{
  double x;
  std::memcpy(&x, &u, sizeof(x));
  return x;
}

uint64_t to_bits(double x)
{
  uint64_t u;
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

double call(const std::string& name, double x)
{
  if (name == "exp")
    return openmc::exp(x);
  if (name == "expm1")
    return openmc::expm1(x);
  if (name == "log")
    return openmc::log(x);
  if (name == "log1p")
    return openmc::log1p(x);
  if (name == "log2")
    return openmc::log2(x);
  if (name == "sin")
    return openmc::sin(x);
  if (name == "cos")
    return openmc::cos(x);
  if (name == "tan")
    return openmc::tan(x);
  if (name == "asin")
    return openmc::asin(x);
  if (name == "acos")
    return openmc::acos(x);
  if (name == "atan")
    return openmc::atan(x);
  if (name == "sinh")
    return openmc::sinh(x);
  if (name == "cosh")
    return openmc::cosh(x);
  if (name == "erf")
    return openmc::erf(x);
  if (name == "erfc")
    return openmc::erfc(x);
  if (name == "lgamma")
    return openmc::lgamma(x);
  if (name == "tgamma")
    return openmc::tgamma(x);
  if (name == "cbrt")
    return openmc::cbrt(x);
  FAIL("unknown unary function " << name);
  return 0.0;
}

double call(const std::string& name, double x, double y)
{
  if (name == "pow")
    return openmc::pow(x, y);
  if (name == "atan2")
    return openmc::atan2(x, y);
  if (name == "hypot")
    return openmc::hypot(x, y);
  FAIL("unknown binary function " << name);
  return 0.0;
}

// NaN payloads and the sign of a default NaN differ between architectures,
// so NaN results are compared as NaN rather than bit for bit.
void check_bits(const std::string& what, double actual, uint64_t expected)
{
  INFO(what);
  if (std::isnan(from_bits(expected))) {
    CHECK(std::isnan(actual));
  } else {
    CHECK(to_bits(actual) == expected);
  }
}

} // namespace

TEST_CASE("CORE-MATH functions reproduce the known-answer vectors")
{
  if (!OPENMC_PORTABLE_MATH_AVAILABLE)
    SKIP("built without OPENMC_USE_COREMATH");

  for (const auto& v : portable_math_kat::unary) {
    double x = from_bits(v.x);
    check_bits(std::string(v.name) + "(" + std::to_string(x) + ")",
      call(v.name, x), v.result);
  }
  for (const auto& v : portable_math_kat::binary) {
    double x = from_bits(v.x);
    double y = from_bits(v.y);
    check_bits(std::string(v.name) + "(" + std::to_string(x) + ", " +
                 std::to_string(y) + ")",
      call(v.name, x, y), v.result);
  }
}

TEST_CASE("CORE-MATH functions handle special values")
{
  if (!OPENMC_PORTABLE_MATH_AVAILABLE)
    SKIP("built without OPENMC_USE_COREMATH");

  const double inf = INFINITY;
  const double nan = NAN;

  CHECK(openmc::exp(0.0) == 1.0);
  CHECK(openmc::exp(-inf) == 0.0);
  CHECK(openmc::exp(inf) == inf);
  CHECK(openmc::exp(710.0) == inf);
  CHECK(openmc::exp(-746.0) == 0.0);
  CHECK(std::isnan(openmc::exp(nan)));

  CHECK(openmc::log(1.0) == 0.0);
  CHECK(openmc::log(0.0) == -inf);
  CHECK(std::isnan(openmc::log(-1.0)));
  CHECK(openmc::log(inf) == inf);
  CHECK(openmc::log2(8.0) == 3.0);
  CHECK(openmc::log1p(0.0) == 0.0);
  CHECK(openmc::expm1(0.0) == 0.0);

  CHECK(openmc::sin(0.0) == 0.0);
  CHECK(std::signbit(openmc::sin(-0.0)));
  CHECK(openmc::cos(0.0) == 1.0);
  CHECK(std::isnan(openmc::sin(inf)));
  CHECK(std::isnan(openmc::cos(inf)));
  CHECK(openmc::tan(0.0) == 0.0);
  CHECK(openmc::asin(1.0) == 0x1.921fb54442d18p+0);
  CHECK(openmc::acos(1.0) == 0.0);
  CHECK(openmc::acos(-1.0) == 0x1.921fb54442d18p+1);
  CHECK(std::isnan(openmc::asin(1.5)));
  CHECK(openmc::atan(inf) == 0x1.921fb54442d18p+0);
  CHECK(openmc::atan2(0.0, -1.0) == 0x1.921fb54442d18p+1);
  CHECK(openmc::atan2(1.0, 0.0) == 0x1.921fb54442d18p+0);

  CHECK(openmc::sinh(0.0) == 0.0);
  CHECK(openmc::cosh(0.0) == 1.0);
  CHECK(openmc::cosh(710.0) == 0x1.3e21a464507f9p+1023);
  CHECK(openmc::cosh(711.0) == inf);
  CHECK(openmc::cosh(-5.0) == openmc::cosh(5.0));

  CHECK(openmc::erf(0.0) == 0.0);
  CHECK(openmc::erf(inf) == 1.0);
  CHECK(openmc::erf(-inf) == -1.0);
  CHECK(openmc::erfc(0.0) == 1.0);
  CHECK(openmc::erfc(inf) == 0.0);
  CHECK(openmc::tgamma(1.0) == 1.0);
  CHECK(openmc::tgamma(5.0) == 24.0);
  CHECK(openmc::lgamma(1.0) == 0.0);
  CHECK(openmc::lgamma(2.0) == 0.0);

  // pow near the overflow and underflow thresholds and with extreme exponents
  CHECK(openmc::pow(1.0, 1e308) == 1.0);
  CHECK(openmc::pow(1.0, nan) == 1.0);
  CHECK(openmc::pow(nan, 0.0) == 1.0);
  CHECK(openmc::pow(2.0, 1e-300) == 1.0);
  CHECK(openmc::pow(0.5, 1e308) == 0.0);
  CHECK(openmc::pow(2.0, 1e308) == inf);
  CHECK(openmc::pow(2.0, -1e308) == 0.0);
  CHECK(openmc::pow(1e300, -1.0) == 1e-300);
  CHECK(openmc::pow(1.7976931348623157e308, -1.0) == 5.562684646268003e-309);
  CHECK(
    openmc::pow(2.718281828459045, 709.782712893384) == 1.7976931348622053e308);
  CHECK(openmc::pow(2.718281828459045, 709.7827128933841) == inf);
  CHECK(openmc::pow(2.718281828459045, -745.1332191019411) == 5e-324);
  CHECK(openmc::pow(0.0, -1.0) == inf);
  CHECK(openmc::pow(-0.0, -1.0) == -inf);
  CHECK(openmc::pow(-8.0, 1.0 / 3.0) != openmc::pow(-8.0, 1.0 / 3.0)); // NaN
  CHECK(openmc::pow(-2.0, 3.0) == -8.0);
  CHECK(openmc::pow(-2.0, 2.0) == 4.0);
  CHECK(openmc::pow(4.0, 0.5) == 2.0);
  CHECK(openmc::pow(3.0, 2.0) == 9.0);

  CHECK(openmc::hypot(3.0, 4.0) == 5.0);
  CHECK(openmc::hypot(inf, nan) == inf);
  CHECK(std::isnan(openmc::hypot(1.0, nan)));
  CHECK(openmc::hypot(0.0, 0.0) == 0.0);
  CHECK(openmc::hypot(1e308, 1e308) == 0x1.92c80954c51f5p+1023);

  CHECK(openmc::cbrt(0.0) == 0.0);
  CHECK(std::signbit(openmc::cbrt(-0.0)));
  CHECK(openmc::cbrt(-8.0) == -2.0);
  CHECK(openmc::cbrt(27.0) == 3.0);
  CHECK(openmc::cbrt(inf) == inf);
  CHECK(openmc::cbrt(-inf) == -inf);
  CHECK(openmc::cbrt(5e-324) == 0x1p-358);
  CHECK(openmc::cbrt(1.7976931348623157e308) == 0x1.428a2f98d728bp+341);
}
