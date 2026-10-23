#include "openmc/scattdata.h"

#include <cmath>
#include <type_traits>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/random_lcg.h"

namespace {

using Catch::Matchers::WithinAbs;
using openmc::double_2dvec;
using openmc::double_3dvec;
using openmc::MgxsType;

// Use the same isotropic production cross sections in all representations.
// With two histogram bins, each bin has width one.
template<typename T>
T make_scatter(const double_2dvec& production, const double_2dvec& multiplicity)
{
  openmc::tensor::Tensor<int> gmin({2}, 0);
  openmc::tensor::Tensor<int> gmax({2}, 1);
  gmin[1] = 1;

  double_3dvec coefficients(2);
  for (int gin = 0; gin < 2; ++gin) {
    for (double value : production[gin]) {
      if constexpr (std::is_same_v<T, openmc::ScattDataLegendre>) {
        coefficients[gin].push_back({value});
      } else {
        coefficients[gin].push_back({value / 2.0, value / 2.0});
      }
    }
  }

  T data;
  data.init(gmin, gmax, multiplicity, coefficients);
  return data;
}

void check_cross_sections(openmc::ScattData& data,
  const double_2dvec& production, const double_2dvec& scattering)
{
  for (int gin = 0; gin < 2; ++gin) {
    CAPTURE(gin);
    double sum_production = 0.0;
    double sum_scattering = 0.0;
    for (int i = 0; i < production[gin].size(); ++i) {
      sum_production += production[gin][i];
      sum_scattering += scattering[gin][i];
      int gout = data.gmin[gin] + i;
      CHECK_THAT(data.get_xs(MgxsType::NU_SCATTER, gin, &gout, nullptr),
        WithinAbs(production[gin][i], 1e-12));
      CHECK_THAT(data.get_xs(MgxsType::SCATTER, gin, &gout, nullptr),
        WithinAbs(scattering[gin][i], 1e-12));
    }
    CHECK_THAT(data.get_xs(MgxsType::NU_SCATTER, gin, nullptr, nullptr),
      WithinAbs(sum_production, 1e-12));
    CHECK_THAT(data.get_xs(MgxsType::SCATTER, gin, nullptr, nullptr),
      WithinAbs(sum_scattering, 1e-12));
  }
}

void check_sample_weights(
  openmc::ScattData& data, int gin, double expected_multiplicity)
{
  uint64_t seed = openmc::init_seed(1, 0);
  constexpr double initial_weight = 2.5;
  for (int i = 0; i < 32; ++i) {
    CAPTURE(gin, i);
    int gout;
    double mu;
    double weight = initial_weight;
    data.sample(gin, gout, mu, weight, &seed);
    REQUIRE(gout >= data.gmin[gin]);
    REQUIRE(gout <= data.gmax[gin]);
    REQUIRE(mu >= -1.0);
    REQUIRE(mu <= 1.0);
    // The weight multiplier is fixed by the incoming group, regardless of the
    // outgoing group drawn from the production-weighted distribution.
    REQUIRE_THAT(
      weight, WithinAbs(initial_weight * expected_multiplicity, 1e-12));
  }
}

} // namespace

TEMPLATE_TEST_CASE(
  "Multigroup scattering cross sections with variable multiplicity",
  "[scattdata]", openmc::ScattDataLegendre, openmc::ScattDataHistogram,
  openmc::ScattDataTabular)
{
  auto data = make_scatter<TestType>({{0.2, 0.4}, {0.2}}, {{1.0, 2.0}, {1.0}});
  CHECK_THAT(data.energy[0][0], WithinAbs(1.0 / 3.0, 1e-12));
  CHECK_THAT(data.energy[0][1], WithinAbs(2.0 / 3.0, 1e-12));
  CHECK_THAT(data.energy[1][0], WithinAbs(1.0, 1e-12));
  check_cross_sections(data, {{0.2, 0.4}, {0.2}}, {{0.2, 0.2}, {0.2}});
}

TEMPLATE_TEST_CASE("Multigroup sampling with variable multiplicity",
  "[scattdata]", openmc::ScattDataLegendre, openmc::ScattDataHistogram,
  openmc::ScattDataTabular)
{
  auto data = make_scatter<TestType>({{0.2, 0.4}, {0.2}}, {{1.0, 2.0}, {1.0}});
  check_sample_weights(data, 0, 1.5);
  check_sample_weights(data, 1, 1.0);
}

TEMPLATE_TEST_CASE("Multigroup sampling with constant multiplicity",
  "[scattdata]", openmc::ScattDataLegendre, openmc::ScattDataHistogram,
  openmc::ScattDataTabular)
{
  double multiplicity = 1.0;
  SECTION("Unit multiplicity") {}
  SECTION("Nonunit multiplicity")
  {
    multiplicity = 1.5;
  }

  auto data = make_scatter<TestType>(
    {{0.2, 0.4}, {0.2}}, {{multiplicity, multiplicity}, {multiplicity}});
  CHECK_THAT(data.energy[0][0], WithinAbs(1.0 / 3.0, 1e-12));
  CHECK_THAT(data.energy[0][1], WithinAbs(2.0 / 3.0, 1e-12));
  check_cross_sections(data, {{0.2, 0.4}, {0.2}},
    {{0.2 / multiplicity, 0.4 / multiplicity}, {0.2 / multiplicity}});
  check_sample_weights(data, 0, multiplicity);
  check_sample_weights(data, 1, multiplicity);
}

TEMPLATE_TEST_CASE("Combining multigroup data with different multiplicities",
  "[scattdata]", openmc::ScattDataLegendre, openmc::ScattDataHistogram,
  openmc::ScattDataTabular)
{
  auto first = make_scatter<TestType>({{0.2, 0.4}, {0.2}}, {{1.0, 2.0}, {1.0}});
  auto second =
    make_scatter<TestType>({{0.3, 0.1}, {0.3}}, {{3.0, 1.0}, {1.5}});
  TestType combined;
  combined.combine({&first, &second}, {0.5, 2.0});

  // Cross sections, rather than multiplicities, add with number density.
  check_cross_sections(combined, {{0.7, 0.4}, {0.7}}, {{0.3, 0.3}, {0.5}});
  CHECK_THAT(combined.energy[0][0], WithinAbs(7.0 / 11.0, 1e-12));
  CHECK_THAT(combined.energy[0][1], WithinAbs(4.0 / 11.0, 1e-12));
  check_sample_weights(combined, 0, 11.0 / 6.0);
  check_sample_weights(combined, 1, 1.4);
}

TEST_CASE("Legendre to tabular conversion preserves variable multiplicity",
  "[scattdata]")
{
  openmc::tensor::Tensor<int> gmin({2}, 0);
  openmc::tensor::Tensor<int> gmax({2}, 1);
  gmin[1] = 1;
  openmc::ScattDataLegendre legendre;
  // A positive first-order angular distribution exercises the conversion
  // without requiring the negative-probability correction.
  legendre.init(gmin, gmax, {{1.0, 2.0}, {1.0}},
    {{{0.2, 0.04}, {0.4, 0.08}}, {{0.2, 0.04}}});

  openmc::ScattDataTabular tabular;
  convert_legendre_to_tabular(legendre, tabular);
  check_cross_sections(tabular, {{0.2, 0.4}, {0.2}}, {{0.2, 0.2}, {0.2}});
  CHECK_THAT(tabular.energy[0][0], WithinAbs(1.0 / 3.0, 1e-12));
  CHECK_THAT(tabular.energy[0][1], WithinAbs(2.0 / 3.0, 1e-12));
  for (double mu : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
    CHECK_THAT(
      tabular.calc_f(0, 0, mu), WithinAbs(legendre.calc_f(0, 0, mu), 1e-12));
    CHECK_THAT(
      tabular.calc_f(0, 1, mu), WithinAbs(legendre.calc_f(0, 1, mu), 1e-12));
  }
  check_sample_weights(tabular, 0, 1.5);
  check_sample_weights(tabular, 1, 1.0);
}

TEMPLATE_TEST_CASE("Zero multigroup scattering row has finite cross sections",
  "[scattdata]", openmc::ScattDataLegendre, openmc::ScattDataHistogram,
  openmc::ScattDataTabular)
{
  auto data = make_scatter<TestType>({{0.0, 0.0}, {0.2}}, {{0.0, 0.0}, {1.0}});
  check_cross_sections(data, {{0.0, 0.0}, {0.2}}, {{0.0, 0.0}, {0.2}});
  CHECK(std::isfinite(data.get_xs(MgxsType::SCATTER, 0, nullptr, nullptr)));
  // This incoming group has no scattering events and must never be sampled.
  check_sample_weights(data, 1, 1.0);
}
