#include <cmath>
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

#include "openmc/math_functions.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/wmp.h"

TEST_CASE("Test t_percentile")
{
  // The reference solutions come from scipy.stats.t.ppf
  std::vector<std::vector<double>> ref_ts {
    {-15.894544844102773, -0.32491969623407446, 0.000000000000000,
      0.32491969623407446, 15.894544844102759},
    {-4.848732214442601, -0.2886751346880066, 0.000000000000000,
      0.2886751346880066, 4.848732214442598},
    {-2.756508521909475, -0.2671808657039658, 0.000000000000000,
      0.2671808657039658, 2.7565085219094745}};

  // Permutations include 1 DoF, 2 DoF, and > 2 DoF
  // We will test 5 p-values at 3-DoF values
  std::vector<double> test_ps {0.02, 0.4, 0.5, 0.6, 0.98};
  std::vector<int> test_dfs {1, 2, 5};

  for (int i = 0; i < test_dfs.size(); i++) {
    int df = test_dfs[i];

    std::vector<double> test_ts;

    for (double p : test_ps) {
      double test_t = openmc::t_percentile(p, df);
      test_ts.push_back(test_t);
    }

    // The 5 DoF approximation in openmc.lib.math.t_percentile is off by up to
    // 8e-3 from the scipy solution, so test that one separately with looser
    // tolerance
    double tolerance = (df > 2) ? 1e-2 : 1e-6;

    REQUIRE_THAT(
      ref_ts[i], Catch::Matchers::Approx(test_ts).epsilon(tolerance));
  }
}

TEST_CASE("Test calc_pn")
{
  // The reference solutions come from scipy.special.eval_legendre
  std::vector<std::vector<double>> ref_vals {
    {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
    {1, -0.5, -0.125, 0.4375, -0.289062, -0.0898438, 0.323242, -0.223145,
      -0.0736389, 0.267899, -0.188229},
    {1, 0, -0.5, -0, 0.375, 0, -0.3125, -0, 0.273438, 0, -0.246094},
    {1, 0.5, -0.125, -0.4375, -0.289062, 0.0898438, 0.323242, 0.223145,
      -0.0736389, -0.267899, -0.188229},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

  int max_order = 10;
  std::vector<double> test_xs = {-1.0, -0.5, 0.0, 0.5, 1.0};

  std::vector<std::vector<double>> test_vals;
  for (double x : test_xs) {
    std::vector<double> test_val(max_order + 1);
    openmc::calc_pn_c(max_order, x, test_val.data());
    test_vals.push_back(test_val);
  }

  for (int i = 0; i < ref_vals.size(); i++) {
    REQUIRE_THAT(ref_vals[i], Catch::Matchers::Approx(test_vals[i]));
  }
}

TEST_CASE("Test evaluate_legendre")
{
  // The reference solutions come from numpy.polynomial.legendre.legval
  std::vector<double> ref_vals {
    5.5, -0.45597649, -1.35351562, -2.7730999, 60.5};

  int max_order = 10;
  std::vector<double> test_xs = {-1.0, -0.5, 0.0, 0.5, 1.0};

  // Set the coefficients back to 1s for the test values since
  // evaluate legendre incorporates the (2l+1)/2 term on its own
  std::vector<double> test_coeffs(max_order + 1, 1.0);

  std::vector<double> test_vals;
  for (double x : test_xs) {
    test_vals.push_back(
      openmc::evaluate_legendre(test_coeffs.size() - 1, test_coeffs.data(), x));
  }

  REQUIRE_THAT(ref_vals, Catch::Matchers::Approx(test_vals));
}

TEST_CASE("Test calc_rn")
{
  std::vector<double> ref_vals {1.000000000000000, -0.019833838076210,
    0.980066577841242, -0.197676811654084, 0.006790834062088,
    -0.033668438114859, 0.940795745502164, -0.335561350977312,
    0.033500236162691, -0.001831975566765, 0.014882082223994,
    -0.046185860057145, 0.883359726009014, -0.460318044571973,
    0.073415616482180, -0.005922278973373, 0.000448625292461,
    -0.004750335422039, 0.025089695062177, -0.057224052171859,
    0.809468042300133, -0.570331780454957, 0.123771351522967,
    -0.015356543011155, 0.001061098599927, -0.000104097571795,
    0.001319047965347, -0.009263463267120, 0.037043163155191,
    -0.066518621473934, 0.721310852552881, -0.662967447756079,
    0.182739660926192, -0.029946258412359, 0.003119841820746,
    -0.000190549327031, 0.000023320052630, -0.000338370521658,
    0.002878809439524, -0.015562587450914, 0.050271226423217,
    -0.073829294593737, 0.621486505922182, -0.735830327235834,
    0.247995745731425, -0.050309614442385, 0.006809024629381,
    -0.000619383085285, 0.000034086826414, -0.000005093626712,
    0.000082405610567, -0.000809532012556, 0.005358034016708,
    -0.023740240859138, 0.064242405926477, -0.078969918083157,
    0.512915839160049, -0.787065093668736, 0.316917738015632,
    -0.076745744765114, 0.012672942183651, -0.001481838409317,
    0.000120451946983, -0.000006047366709, 0.000001091052697,
    -0.000019334294214, 0.000213051604838, -0.001640234119608,
    0.008982263900105, -0.033788039035668, 0.078388909900756,
    -0.081820779415058, 0.398746190829636, -0.815478614863816,
    0.386704633068855, -0.109227544713261, 0.021245051959237,
    -0.003002428416676, 0.000311416667310, -0.000022954482885,
    0.000001059646310, -0.000000230023931, 0.000004408854505,
    -0.000053457925526, 0.000464152759861, -0.002976305522860,
    0.013958017448970, -0.045594791382625, 0.092128969315914,
    -0.082334538374971, 0.282248459574595, -0.820599067736528,
    0.454486474163594, -0.147395565311743, 0.033013815809602,
    -0.005448090715661, 0.000678450207914, -0.000063467485444,
    0.000004281943868, -0.000000182535754, 0.000000047847775,
    -0.000000982664801, 0.000012933320414, -0.000124076425457,
    0.000901739739837, -0.004982323311961, 0.020457776068931,
    -0.058948376674391, 0.104888993733747, -0.080538298991650,
    0.166710763818175, -0.802696588503912, 0.517433650833039,
    -0.190564076304612, 0.048387190622376, -0.009120081648146,
    0.001318069323039, -0.000147308722683, 0.000012561029621,
    -0.000000779794781, 0.000000030722703};

  int max_order = 10;

  double azi = 0.1; // Longitude
  double pol = 0.2; // Latitude
  double mu = std::cos(pol);

  std::vector<double> test_uvw {std::sin(pol) * std::cos(azi),
    std::sin(pol) * std::sin(azi), std::cos(pol)};

  std::vector<double> test_vals((max_order + 1) * (max_order + 1), 0);
  openmc::calc_rn_c(max_order, test_uvw.data(), test_vals.data());

  REQUIRE_THAT(ref_vals, Catch::Matchers::Approx(test_vals));
}

TEST_CASE("Test calc_zn")
{
  std::vector<double> ref_vals {1.00000000e+00, 2.39712769e-01, 4.38791281e-01,
    2.10367746e-01, -5.00000000e-01, 1.35075576e-01, 1.24686873e-01,
    -2.99640962e-01, -5.48489101e-01, 8.84215021e-03, 5.68310892e-02,
    -4.20735492e-01, -1.25000000e-01, -2.70151153e-01, -2.60091773e-02,
    1.87022545e-02, -3.42888902e-01, 1.49820481e-01, 2.74244551e-01,
    -2.43159131e-02, -2.50357380e-02, 2.20500013e-03, -1.98908812e-01,
    4.07587508e-01, 4.37500000e-01, 2.61708929e-01, 9.10321205e-02,
    -1.54686328e-02, -2.74049397e-03, -7.94845816e-02, 4.75368705e-01,
    7.11647284e-02, 1.30266162e-01, 3.37106977e-02, 1.06401886e-01,
    -7.31606787e-03, -2.95625975e-03, -1.10250006e-02, 3.55194307e-01,
    -1.44627826e-01, -2.89062500e-01, -9.28644588e-02, -1.62557358e-01,
    7.73431638e-02, -2.55329539e-03, -1.90923851e-03, 1.57578403e-02,
    1.72995854e-01, -3.66267690e-01, -1.81657333e-01, -3.32521518e-01,
    -2.59738162e-02, -2.31580576e-01, 4.20673902e-02, -4.11710546e-04,
    -9.36449487e-04, 1.92156884e-02, 2.82515641e-02, -3.90713738e-01,
    -1.69280296e-01, -8.98437500e-02, -1.08693628e-01, 1.78813094e-01,
    -1.98191857e-01, 1.65964201e-02, 2.77013853e-04};

  int n = 10;
  double rho = 0.5;
  double phi = 0.5;

  int nums = ((n + 1) * (n + 2)) / 2;

  std::vector<double> test_vals(nums, 0);
  openmc::calc_zn(n, rho, phi, test_vals.data());

  REQUIRE_THAT(ref_vals, Catch::Matchers::Approx(test_vals));
}

TEST_CASE("Test calc_zn_rad")
{
  std::vector<double> ref_vals {1.00000000e+00, -5.00000000e-01,
    -1.25000000e-01, 4.37500000e-01, -2.89062500e-01, -8.98437500e-02};

  int n = 10;
  double rho = 0.5;

  int nums = n / 2 + 1;
  std::vector<double> test_vals(nums, 0);
  openmc::calc_zn_rad(n, rho, test_vals.data());

  REQUIRE_THAT(ref_vals, Catch::Matchers::Approx(test_vals));
}

TEST_CASE("Test rotate_angle")
{
  std::vector<double> uvw0 {1.0, 0.0, 0.0};
  double phi = 0.0;

  uint64_t prn_seed = 1;
  openmc::prn(&prn_seed);

  SECTION("Test rotate_angle mu is 0")
  {
    std::vector<double> ref_uvw {0.0, 0.0, -1.0};

    double mu = 0.0;

    std::vector<double> test_uvw(uvw0);
    openmc::rotate_angle_c(test_uvw.data(), mu, &phi, &prn_seed);

    REQUIRE_THAT(ref_uvw, Catch::Matchers::Approx(test_uvw));
  }

  SECTION("Test rotate_angle mu is 1")
  {
    std::vector<double> ref_uvw = {1.0, 0.0, 0.0};

    double mu = 1.0;

    std::vector<double> test_uvw(uvw0);
    openmc::rotate_angle_c(test_uvw.data(), mu, &phi, &prn_seed);

    REQUIRE_THAT(ref_uvw, Catch::Matchers::Approx(test_uvw));
  }

  // Now to test phi is None
  SECTION("Test rotate_angle no phi")
  {
    // When seed = 1, phi will be sampled as 1.9116495709698769
    // The resultant reference is from hand-calculations given the above
    std::vector<double> ref_uvw = {
      0.9, -0.422746750548505, 0.10623175090659095};

    double mu = 0.9;
    prn_seed = 1;

    std::vector<double> test_uvw(uvw0);
    openmc::rotate_angle_c(test_uvw.data(), mu, NULL, &prn_seed);

    REQUIRE_THAT(ref_uvw, Catch::Matchers::Approx(test_uvw));
  }
}

TEST_CASE("Test maxwell_spectrum")
{
  double ref_val = 0.27767406743161277;

  double T = 0.5;
  uint64_t prn_seed = 1;

  double test_val = openmc::maxwell_spectrum(T, &prn_seed);

  REQUIRE(ref_val == test_val);
}

TEST_CASE("Test watt_spectrum")
{
  double ref_val = 0.30957476387766697;

  double a = 0.5;
  double b = 0.75;
  uint64_t prn_seed = 1;

  double test_val = openmc::watt_spectrum(a, b, &prn_seed);

  REQUIRE(ref_val == test_val);
}

TEST_CASE("Test normal_variate")
{

  // Generate a series of normally distributed random numbers and test
  // whether their mean and standard deviation are close to the expected value
  SECTION("Test with non-zero standard deviation")
  {
    uint64_t seed = 1;

    double mean = 0.0;
    double standard_deviation = 1.0;

    int num_samples = 10000;
    double sum = 0.0;
    double sum_squared_difference = 0.0;

    for (int i = 0; i < num_samples; ++i) {
      double sample = openmc::normal_variate(mean, standard_deviation, &seed);
      sum += sample;
      sum_squared_difference += (sample - mean) * (sample - mean);
    }

    double actual_mean = sum / num_samples;
    double actual_standard_deviation =
      std::sqrt(sum_squared_difference / num_samples);

    REQUIRE_THAT(mean, Catch::Matchers::WithinAbs(actual_mean, 0.1));
    REQUIRE_THAT(standard_deviation,
      Catch::Matchers::WithinAbs(actual_standard_deviation, 0.1));
  }

  // When the standard deviation is zero
  // the generated random number should always be equal to the mean
  SECTION("Test with zero standard deviation")
  {
    uint64_t seed = 1;
    double mean = 5.0;
    double standard_deviation = 0.0;

    for (int i = 0; i < 10; ++i) {
      double sample = openmc::normal_variate(mean, standard_deviation, &seed);
      REQUIRE(sample == mean);
    }
  }
}

TEST_CASE("Test broaden_wmp_polynomials")
{
  double test_E = 0.5;
  int n = 6;

  // Two branches of the code to worry about, beta > 6 and otherwise
  // beta = sqrtE * dopp
  SECTION("Test broaden_wmp_polynomials beta > 6")
  {
    std::vector<double> ref_val {
      2., 1.41421356, 1.0001, 0.70731891, 0.50030001, 0.353907};

    double test_dopp = 100.0; // approximately U235 at room temperature

    std::vector<double> test_val(n, 0);
    openmc::broaden_wmp_polynomials(test_E, test_dopp, n, test_val.data());

    REQUIRE_THAT(ref_val, Catch::Matchers::Approx(test_val));
  }

  SECTION("Test broaden_wmp_polynomials beta < 6")
  {
    std::vector<double> ref_val = {
      1.99999885, 1.41421356, 1.04, 0.79195959, 0.6224, 0.50346003};

    double test_dopp = 5.0;

    std::vector<double> test_val(n, 0);
    openmc::broaden_wmp_polynomials(test_E, test_dopp, n, test_val.data());

    REQUIRE_THAT(ref_val, Catch::Matchers::Approx(test_val));
  }
}
