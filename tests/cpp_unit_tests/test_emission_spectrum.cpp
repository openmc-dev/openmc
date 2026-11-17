#include "openmc/capi.h"
#include "openmc/distribution_energy.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/simulation.h"

#include <cmath>
#include <limits>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

namespace {
using openmc::vector;

// Exercise the real HDF5 loader without an external nuclear-data library.
openmc::ContinuousTabular spectrum(int incident_interpolation = 2)
{
  hid_t properties = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_core(properties, 4096, false);
  hid_t file =
    H5Fcreate("emission-test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, properties);
  H5Pclose(properties);
  openmc::write_dataset(file, "energy", vector<double> {1.0, 3.0});
  hid_t energy = openmc::open_dataset(file, "energy");
  int interp[] {2, incident_interpolation};
  hsize_t shape[] {2, 1};
  openmc::write_attr(energy, 2, shape, "interpolation", H5T_NATIVE_INT, interp);
  openmc::close_dataset(energy);
  auto values = openmc::tensor::zeros<double>({3, 4});
  for (int i = 0; i < 2; ++i) {
    values(0, 2 * i) = 0.0;
    values(0, 2 * i + 1) = 2.0;
    values(1, 2 * i) = values(1, 2 * i + 1) = 0.5;
    values(2, 2 * i + 1) = 1.0;
  }
  openmc::write_dataset(file, "distribution", values);
  hid_t table = openmc::open_dataset(file, "distribution");
  openmc::write_attribute(table, "offsets", vector<int> {0, 2});
  openmc::write_attribute(table, "interpolation", vector<int> {2, 2});
  openmc::write_attribute(table, "n_discrete_lines", vector<int> {0, 0});
  openmc::close_dataset(table);
  openmc::ContinuousTabular result(file);
  H5Fclose(file);
  return result;
}

struct SimulationState {
  bool initialized = openmc::simulation::initialized;
  bool in_batch = openmc::simulation::in_batch;
  ~SimulationState()
  {
    openmc::simulation::initialized = initialized;
    openmc::simulation::in_batch = in_batch;
  }
};
} // namespace

TEST_CASE("Continuous-tabular replacement normalizes and samples mixed tables")
{
  int interpolation = GENERATE(1, 2);
  auto dist = spectrum();
  REQUIRE(dist.energy() == vector<double> {1.0, 3.0});
  REQUIRE(dist.distribution()[0].p[0] == 0.5);
  dist.set_spectrum({1.0, 3.0}, {interpolation, interpolation}, {1, 1},
    {{4.0, 0.0, 2.0}, {4.0, 0.0, 2.0}}, {{2.0, 1.0, 1.0}, {2.0, 1.0, 1.0}});
  REQUIRE(dist.distribution()[0].c[0] == 0.5);
  REQUIRE(dist.distribution()[0].c[1] == 0.5);
  REQUIRE(dist.distribution()[0].c[2] == 1.0);
  uint64_t seed = openmc::init_seed(1, 0);
  double sum = 0.0;
  int lines = 0;
  for (int i = 0; i < 20000; ++i) {
    double value = dist.sample(2.0, &seed);
    REQUIRE((value == 4.0 || (value >= 0.0 && value <= 2.0)));
    sum += value;
    lines += value == 4.0;
  }
  REQUIRE(std::abs(sum / 20000 - 2.5) < 0.05);
  REQUIRE(std::abs(lines / 20000.0 - 0.5) < 0.02);
}

TEST_CASE("Discrete-only tables never access a continuous interval")
{
  int interpolation = GENERATE(1, 2);
  auto dist = spectrum();
  dist.set_spectrum({1.0, 3.0}, {interpolation, interpolation}, {1, 2},
    {{4.0}, {2.0, 7.0}}, {{3.0}, {0.0, 2.0}});
  uint64_t seed = openmc::init_seed(2, 0);
  for (int i = 0; i < 1000; ++i) {
    REQUIRE(dist.sample(0.5, &seed) == 4.0);
    REQUIRE(dist.sample(4.0, &seed) == 7.0);
    double value = dist.sample(2.0, &seed);
    REQUIRE((value == 4.0 || value == 7.0));
  }
}

TEST_CASE("Continuous tables sample zero-probability intervals safely")
{
  int interpolation = GENERATE(1, 2);
  auto dist = spectrum();
  dist.set_spectrum({1.0, 3.0}, {interpolation, interpolation}, {0, 0},
    {{0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 2.0, 3.0}},
    {{0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 1.0, 0.0}});
  uint64_t seed = openmc::init_seed(3, 0);
  double sum = 0.0;
  for (int i = 0; i < 20000; ++i) {
    double value = dist.sample(2.0, &seed);
    REQUIRE(std::isfinite(value));
    REQUIRE(value >= (interpolation == 1 ? 2.0 : 1.0));
    REQUIRE(value <= 3.0);
    sum += value;
  }
  REQUIRE(std::abs(sum / 20000 - (interpolation == 1 ? 2.5 : 2.0)) < 0.03);
}

TEST_CASE("Invalid spectrum replacements leave every native table unchanged")
{
  auto dist = spectrum();
  vector<double> grid {1.0, 3.0};
  vector<int> interp {2, 2}, discrete {0, 0};
  vector<vector<double>> energy {{0.0, 2.0}, {0.0, 2.0}};
  vector<vector<double>> pdf {{0.0, 2.0}, {1.0, 1.0}};
  SECTION("changed grid")
  {
    grid[1] = 4.0;
  }
  SECTION("short grid")
  {
    grid.pop_back();
  }
  SECTION("missing table")
  {
    energy.pop_back();
  }
  SECTION("invalid interpolation")
  {
    interp[1] = 3;
  }
  SECTION("negative discrete count")
  {
    discrete[1] = -1;
  }
  SECTION("excess discrete count")
  {
    discrete[1] = 3;
  }
  SECTION("one continuous point")
  {
    discrete[1] = 1;
  }
  SECTION("negative energy")
  {
    energy[1][0] = -1.0;
  }
  SECTION("nonfinite energy")
  {
    energy[1][0] = INFINITY;
  }
  SECTION("unordered energy")
  {
    energy[1] = {2.0, 0.0};
  }
  SECTION("duplicate energy")
  {
    energy[1] = {2.0, 2.0};
  }
  SECTION("unequal lengths")
  {
    pdf[1].pop_back();
  }
  SECTION("empty table")
  {
    energy[1].clear();
    pdf[1].clear();
  }
  SECTION("negative PDF")
  {
    pdf[1][1] = -1.0;
  }
  SECTION("nonfinite PDF")
  {
    pdf[1][1] = NAN;
  }
  SECTION("zero normalization")
  {
    pdf[1] = {0.0, 0.0};
  }
  SECTION("overflow normalization")
  {
    pdf[1] = {1.e308, 1.e308};
  }
  SECTION("incompatible neighboring table")
  {
    discrete[1] = 2;
  }
  REQUIRE_THROWS_AS(dist.set_spectrum(grid, interp, discrete, energy, pdf),
    std::invalid_argument);
  for (const auto& table : dist.distribution()) {
    REQUIRE(table.e_out.size() == 2);
    REQUIRE(table.e_out[1] == 2.0);
    REQUIRE(table.p[0] == 0.5);
    REQUIRE(table.p[1] == 0.5);
    REQUIRE(table.c[0] == 0.0);
    REQUIRE(table.c[1] == 1.0);
  }
}

TEST_CASE("Histogram incident interpolation permits independent table types")
{
  auto dist = spectrum(1);
  dist.set_spectrum(
    {1.0, 3.0}, {2, 2}, {0, 1}, {{0.0, 2.0}, {7.0}}, {{1.0, 1.0}, {1.0}});
  uint64_t seed = openmc::init_seed(4, 0);
  for (int i = 0; i < 100; ++i) {
    double value = dist.sample(2.0, &seed);
    REQUIRE(value >= 0.0);
    REQUIRE(value <= 2.0);
  }
}

TEST_CASE(
  "Runtime nuclear-data mutation requires an initialized idle simulation")
{
  SimulationState restore;
  bool initialized = GENERATE(false, true);
  openmc::simulation::initialized = initialized;
  openmc::simulation::in_batch = initialized;
  int expected = initialized ? OPENMC_E_INVALID_ARGUMENT : OPENMC_E_ALLOCATE;
  REQUIRE(openmc_nuclide_set_emission_spectrum(0, 18, 0, 0, nullptr, 0, nullptr,
            nullptr, nullptr, nullptr, nullptr, 0) == expected);
  REQUIRE(openmc_nuclide_set_reaction_xs_with_threshold(
            0, 18, 0, nullptr, nullptr, 0, 0.0) == expected);
  REQUIRE(openmc_nuclide_rebuild_derived_xs(0) == expected);
}
