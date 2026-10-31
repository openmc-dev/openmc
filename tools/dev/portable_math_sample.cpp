// Sample the transcendental functions in openmc/math.h (the CORE-MATH backend
// when compiled with -DOPENMC_COREMATH) and the system libm on the same inputs
// and record every result as raw IEEE 754 bits.
//
// The output is consumed by tools/dev/portable_math_report.py, which compares
// both implementations against correctly rounded MPFR references and writes an
// accuracy report. Keeping the sampling in C++ and the reference evaluation in
// Python means the C++ side has no dependency beyond the compiler.
//
// Build (from the repository root):
//
//   g++ -O2 -std=c++17 -DOPENMC_COREMATH -Iinclude
//       tools/dev/portable_math_sample.cpp <build>/lib/libcoremath.a
//       -o portable_math_sample
//
// (one command line; <build> is a build directory configured with
// OPENMC_USE_COREMATH=ON)
//
// Run:
//
//   ./portable_math_sample --points 1000000 --seed 1 > samples.txt
//
// The first line records how many special values follow for each function.
// Each other line is: function, input bits (one or two hex words), openmc::
// result bits, system libm result bits. Inputs are sampled uniformly over the
// IEEE 754 bit patterns of each range, so every representable double in a range
// is equally likely and subnormals, huge arguments and arguments near special
// points are all exercised. A fixed set of special values is emitted first.

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "openmc/math.h"

namespace {

uint64_t to_bits(double x)
{
  uint64_t u;
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

double from_bits(uint64_t u)
{
  double x;
  std::memcpy(&x, &u, sizeof(x));
  return x;
}

// Uniform over the bit patterns of [lo, hi], 0 <= lo <= hi.
double uniform_bits(std::mt19937_64& rng, double lo, double hi)
{
  uint64_t a = to_bits(lo), b = to_bits(hi);
  std::uniform_int_distribution<uint64_t> dist(a, b);
  return from_bits(dist(rng));
}

// Uniform over the bit patterns of [lo, hi] in magnitude, random sign.
double uniform_bits_signed(std::mt19937_64& rng, double lo, double hi)
{
  double x = uniform_bits(rng, lo, hi);
  return (rng() & 1) ? -x : x;
}

struct Range {
  double lo;
  double hi;
  bool random_sign;
  double weight;
};

struct Unary {
  const char* name;
  double (*ours)(double);
  double (*sys)(double);
  std::vector<Range> ranges;
};

struct Binary {
  const char* name;
  double (*ours)(double, double);
  double (*sys)(double, double);
  std::vector<Range> xranges;
  std::vector<Range> yranges;
};

// Wrappers with plain function-pointer types
double o_exp(double x)
{
  return openmc::exp(x);
}
double o_expm1(double x)
{
  return openmc::expm1(x);
}
double o_log(double x)
{
  return openmc::log(x);
}
double o_log1p(double x)
{
  return openmc::log1p(x);
}
double o_log2(double x)
{
  return openmc::log2(x);
}
double o_sin(double x)
{
  return openmc::sin(x);
}
double o_cos(double x)
{
  return openmc::cos(x);
}
double o_tan(double x)
{
  return openmc::tan(x);
}
double o_asin(double x)
{
  return openmc::asin(x);
}
double o_acos(double x)
{
  return openmc::acos(x);
}
double o_atan(double x)
{
  return openmc::atan(x);
}
double o_sinh(double x)
{
  return openmc::sinh(x);
}
double o_cosh(double x)
{
  return openmc::cosh(x);
}
double o_erf(double x)
{
  return openmc::erf(x);
}
double o_erfc(double x)
{
  return openmc::erfc(x);
}
double o_lgamma(double x)
{
  return openmc::lgamma(x);
}
double o_tgamma(double x)
{
  return openmc::tgamma(x);
}
double o_cbrt(double x)
{
  return openmc::cbrt(x);
}
double o_pow(double x, double y)
{
  return openmc::pow(x, y);
}
double o_atan2(double y, double x)
{
  return openmc::atan2(y, x);
}
double o_hypot(double x, double y)
{
  return openmc::hypot(x, y);
}

double s_exp(double x)
{
  return std::exp(x);
}
double s_expm1(double x)
{
  return std::expm1(x);
}
double s_log(double x)
{
  return std::log(x);
}
double s_log1p(double x)
{
  return std::log1p(x);
}
double s_log2(double x)
{
  return std::log2(x);
}
double s_sin(double x)
{
  return std::sin(x);
}
double s_cos(double x)
{
  return std::cos(x);
}
double s_tan(double x)
{
  return std::tan(x);
}
double s_asin(double x)
{
  return std::asin(x);
}
double s_acos(double x)
{
  return std::acos(x);
}
double s_atan(double x)
{
  return std::atan(x);
}
double s_sinh(double x)
{
  return std::sinh(x);
}
double s_cosh(double x)
{
  return std::cosh(x);
}
double s_erf(double x)
{
  return std::erf(x);
}
double s_erfc(double x)
{
  return std::erfc(x);
}
double s_lgamma(double x)
{
  return std::lgamma(x);
}
double s_tgamma(double x)
{
  return std::tgamma(x);
}
double s_cbrt(double x)
{
  return std::cbrt(x);
}
double s_pow(double x, double y)
{
  return std::pow(x, y);
}
double s_atan2(double y, double x)
{
  return std::atan2(y, x);
}
double s_hypot(double x, double y)
{
  return std::hypot(x, y);
}

const double SUBNORMAL_MIN = 4.9406564584124654e-324;
const double TWO_PI = 6.283185307179586;

// Special values every function is evaluated at
const std::vector<double> SPECIALS = {0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 2.0,
  -2.0, 10.0, 0.1, 1e-300, -1e-300, DBL_MIN, -DBL_MIN, SUBNORMAL_MIN,
  -SUBNORMAL_MIN, DBL_MAX, -DBL_MAX, DBL_EPSILON, 1.0 - DBL_EPSILON / 2,
  1.0 + DBL_EPSILON, 0.6931471805599453, 2.718281828459045, 1.5707963267948966,
  3.141592653589793, -3.141592653589793, TWO_PI, 709.782712893384,
  709.7827128933841, -745.1332191019411, -708.3964185322641, 1e22, -1e22, 1e300,
  0.75, 171.62, 171.7, -0.5, -170.5, 0.9999999999999999, INFINITY, -INFINITY,
  NAN};

double pick(std::mt19937_64& rng, const std::vector<Range>& ranges)
{
  double total = 0;
  for (auto& r : ranges)
    total += r.weight;
  std::uniform_real_distribution<double> u(0.0, total);
  double t = u(rng);
  for (auto& r : ranges) {
    if (t < r.weight)
      return r.random_sign ? uniform_bits_signed(rng, r.lo, r.hi)
                           : uniform_bits(rng, r.lo, r.hi);
    t -= r.weight;
  }
  auto& r = ranges.back();
  return r.random_sign ? uniform_bits_signed(rng, r.lo, r.hi)
                       : uniform_bits(rng, r.lo, r.hi);
}

void emit(const Unary& f, double x)
{
  std::printf("%s %016llx %016llx %016llx\n", f.name,
    (unsigned long long)to_bits(x), (unsigned long long)to_bits(f.ours(x)),
    (unsigned long long)to_bits(f.sys(x)));
}

void emit(const Binary& f, double x, double y)
{
  std::printf("%s %016llx %016llx %016llx %016llx\n", f.name,
    (unsigned long long)to_bits(x), (unsigned long long)to_bits(y),
    (unsigned long long)to_bits(f.ours(x, y)),
    (unsigned long long)to_bits(f.sys(x, y)));
}

} // namespace

int main(int argc, char** argv)
{
  size_t points = 100000;
  uint64_t seed = 1;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--points" && i + 1 < argc)
      points = std::strtoull(argv[++i], nullptr, 10);
    else if (a == "--seed" && i + 1 < argc)
      seed = std::strtoull(argv[++i], nullptr, 10);
    else {
      std::fprintf(stderr, "usage: %s [--points N] [--seed S]\n", argv[0]);
      return 2;
    }
  }

  // Whole positive range including subnormals
  const Range all_pos {SUBNORMAL_MIN, DBL_MAX, false, 1.0};
  const Range all_signed {SUBNORMAL_MIN, DBL_MAX, true, 1.0};

  std::vector<Unary> unary = {
    {"exp", o_exp, s_exp,
      {{SUBNORMAL_MIN, 750.0, true, 1.0}, {1e-3, 1.0, true, 1.0}}},
    {"expm1", o_expm1, s_expm1,
      {{SUBNORMAL_MIN, 750.0, true, 1.0}, {1e-3, 1.0, true, 1.0}}},
    {"log", o_log, s_log, {all_pos, {0.5, 2.0, false, 1.0}}},
    {"log1p", o_log1p, s_log1p,
      {{SUBNORMAL_MIN, 1.0, true, 1.0}, {1.0, DBL_MAX, false, 1.0}}},
    {"log2", o_log2, s_log2, {all_pos, {0.5, 2.0, false, 1.0}}},
    {"sin", o_sin, s_sin,
      {{SUBNORMAL_MIN, TWO_PI, true, 2.0}, {TWO_PI, 1e22, true, 1.0}}},
    {"cos", o_cos, s_cos,
      {{SUBNORMAL_MIN, TWO_PI, true, 2.0}, {TWO_PI, 1e22, true, 1.0}}},
    {"tan", o_tan, s_tan,
      {{SUBNORMAL_MIN, TWO_PI, true, 2.0}, {TWO_PI, 1e22, true, 1.0}}},
    {"asin", o_asin, s_asin, {{SUBNORMAL_MIN, 1.0, true, 1.0}}},
    {"acos", o_acos, s_acos, {{SUBNORMAL_MIN, 1.0, true, 1.0}}},
    {"atan", o_atan, s_atan, {all_signed, {0.1, 10.0, true, 1.0}}},
    {"sinh", o_sinh, s_sinh,
      {{SUBNORMAL_MIN, 750.0, true, 1.0}, {1e-3, 1.0, true, 1.0}}},
    {"cosh", o_cosh, s_cosh,
      {{SUBNORMAL_MIN, 750.0, true, 1.0}, {1e-3, 1.0, true, 1.0}}},
    {"erf", o_erf, s_erf,
      {{SUBNORMAL_MIN, 30.0, true, 1.0}, {1e-3, 3.0, true, 1.0}}},
    {"erfc", o_erfc, s_erfc,
      {{SUBNORMAL_MIN, 30.0, true, 1.0}, {1e-3, 3.0, true, 1.0}}},
    {"lgamma", o_lgamma, s_lgamma,
      {{SUBNORMAL_MIN, 1e300, true, 1.0}, {0.5, 200.0, true, 1.0}}},
    {"tgamma", o_tgamma, s_tgamma,
      {{SUBNORMAL_MIN, 180.0, true, 1.0}, {0.5, 50.0, true, 1.0}}},
    {"cbrt", o_cbrt, s_cbrt, {all_signed, {0.1, 1e9, true, 1.0}}},
  };

  // pow inputs: x over the positive range and in the OpenMC energy range
  // (1e-5 eV to 2e7 eV), y over the exponent range where results are finite
  // and the small integer and fractional exponents the code uses.
  std::vector<Binary> binary = {
    {"pow", o_pow, s_pow,
      {all_pos, {1e-5, 2e7, false, 1.0}, {0.5, 2.0, false, 1.0}},
      {{SUBNORMAL_MIN, 1100.0, true, 1.0}, {0.25, 4.0, true, 1.0}}},
    {"atan2", o_atan2, s_atan2, {all_signed, {0.1, 10.0, true, 1.0}},
      {all_signed, {0.1, 10.0, true, 1.0}}},
    {"hypot", o_hypot, s_hypot, {all_signed, {0.1, 1e3, true, 1.0}},
      {all_signed, {0.1, 1e3, true, 1.0}}},
  };

  std::mt19937_64 rng(seed);

  // Lets the report tell special values from sampled ones
  std::printf("# specials %zu\n", SPECIALS.size());

  for (auto& f : unary) {
    for (double x : SPECIALS)
      emit(f, x);
    for (size_t i = 0; i < points; ++i)
      emit(f, pick(rng, f.ranges));
  }
  for (auto& f : binary) {
    for (double x : SPECIALS)
      for (double y : SPECIALS)
        emit(f, x, y);
    // Exponents used in the transport code
    const double fixed_y[] = {2.0, 3.0, -1.0, -2.0, 0.5, 1.0 / 3.0, 1.5};
    for (size_t i = 0; i < points; ++i) {
      double x = pick(rng, f.xranges);
      double y = pick(rng, f.yranges);
      emit(f, x, y);
      if (std::strcmp(f.name, "pow") == 0 && (i % 8) == 0)
        emit(f, x, fixed_y[(i / 8) % 7]);
    }
  }
  return 0;
}
