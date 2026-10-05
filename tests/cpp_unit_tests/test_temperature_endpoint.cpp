#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

namespace {

using namespace openmc;

class H5Handle {
public:
  H5Handle(hid_t id, herr_t (*close)(hid_t)) : id_ {id}, close_ {close}
  {
    if (id < 0)
      throw std::runtime_error("Cannot create temperature endpoint fixture");
  }

  ~H5Handle() { close_(id_); }

  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;

  operator hid_t() const { return id_; }

private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

class TemporaryDirectory {
public:
  TemporaryDirectory()
  {
    std::random_device random;
    for (int attempt = 0; attempt < 16; ++attempt) {
      auto candidate =
        std::filesystem::temp_directory_path() /
        ("openmc-temperature-endpoint-" + std::to_string(random()) + "-" +
          std::to_string(random()));
      if (std::filesystem::create_directory(candidate)) {
        path = std::move(candidate);
        return;
      }
    }
    throw std::runtime_error("Cannot reserve temperature endpoint directory");
  }

  ~TemporaryDirectory()
  {
    std::error_code error;
    std::filesystem::remove_all(path, error);
  }

  TemporaryDirectory(const TemporaryDirectory&) = delete;
  TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

  std::filesystem::path path;
};

struct GlobalState {
  decltype(data::nuclides) nuclides {std::move(data::nuclides)};
  decltype(data::nuclide_map) nuclide_map {std::move(data::nuclide_map)};
  decltype(data::energy_min) energy_min {data::energy_min};
  decltype(data::energy_max) energy_max {data::energy_max};
  double temperature_min {data::temperature_min};
  double temperature_max {data::temperature_max};
  int n_log_bins {settings::n_log_bins};
  RunMode run_mode {settings::run_mode};
  TemperatureMethod temperature_method {settings::temperature_method};
  array<double, 2> temperature_range {settings::temperature_range};
  double temperature_tolerance {settings::temperature_tolerance};
  bool urr_ptables_on {settings::urr_ptables_on};
  bool need_depletion_rx {simulation::need_depletion_rx};

  GlobalState()
  {
    data::nuclides.clear();
    data::nuclide_map.clear();
    data::energy_min = {0.0, 0.0, 0.0, 0.0};
    data::energy_max = {INFTY, INFTY, INFTY, INFTY};
    data::temperature_min = INFTY;
    data::temperature_max = 0.0;
    settings::n_log_bins = 1;
    settings::run_mode = RunMode::FIXED_SOURCE;
    settings::temperature_method = TemperatureMethod::INTERPOLATION;
    settings::temperature_range = {0.0, 0.0};
    settings::temperature_tolerance = 1.0;
    settings::urr_ptables_on = false;
    simulation::need_depletion_rx = false;
  }

  ~GlobalState()
  {
    data::nuclides.clear();
    data::nuclides = std::move(nuclides);
    data::nuclide_map = std::move(nuclide_map);
    data::energy_min = energy_min;
    data::energy_max = energy_max;
    data::temperature_min = temperature_min;
    data::temperature_max = temperature_max;
    settings::n_log_bins = n_log_bins;
    settings::run_mode = run_mode;
    settings::temperature_method = temperature_method;
    settings::temperature_range = temperature_range;
    settings::temperature_tolerance = temperature_tolerance;
    settings::urr_ptables_on = urr_ptables_on;
    simulation::need_depletion_rx = need_depletion_rx;
  }

  GlobalState(const GlobalState&) = delete;
  GlobalState& operator=(const GlobalState&) = delete;
};

int label_for_kT(double kT)
{
  return static_cast<int>(std::lround(kT / K_BOLTZMANN));
}

void write_temperature_table(
  hid_t kts, hid_t energy, hid_t rx, double kT, const vector<double>& xs)
{
  const auto label = std::to_string(label_for_kT(kT)) + "K";
  write_dataset(kts, label.c_str(), kT);
  write_dataset(
    energy, label.c_str(), vector<double> {1.0, 2.0, 3.0, 4.0, 5.0});

  H5Handle temp(create_group(rx, label.c_str()), H5Gclose);
  write_dataset(temp, "xs", xs);
  H5Handle xs_dset(H5Dopen2(temp, "xs", H5P_DEFAULT), H5Dclose);
  write_attribute(xs_dset, "threshold_idx", 0);
}

void load_temperature_fixture(const vector<double>& kTs)
{
  static constexpr double xs_base = 10.0;

  TemporaryDirectory directory;
  const auto pointwise = directory.path / "pointwise.h5";

  H5Handle file(H5Fcreate(pointwise.string().c_str(), H5F_ACC_EXCL, H5P_DEFAULT,
                  H5P_DEFAULT),
    H5Fclose);
  H5Handle group(create_group(file, "test"), H5Gclose);
  write_attribute(group, "Z", 1);
  write_attribute(group, "A", 1);
  write_attribute(group, "metastable", 0);
  write_attribute(group, "atomic_weight_ratio", 1.0);

  H5Handle kts(create_group(group, "kTs"), H5Gclose);
  H5Handle energy(create_group(group, "energy"), H5Gclose);
  H5Handle rxs(create_group(group, "reactions"), H5Gclose);
  H5Handle rx(
    create_group(rxs, "reaction_" + std::to_string(N_GAMMA)), H5Gclose);
  write_attribute(rx, "mt", static_cast<int>(N_GAMMA));
  write_attribute(rx, "Q_value", 0.0);
  write_attribute(rx, "center_of_mass", 0);
  write_attribute(rx, "redundant", 0);

  for (std::size_t i = 0; i < kTs.size(); ++i) {
    const double xs_value = xs_base + static_cast<double>(i);
    write_temperature_table(kts, energy, rx, kTs[i],
      vector<double> {xs_value, xs_value, xs_value, xs_value, xs_value});
  }

  data::nuclides.push_back(make_unique<Nuclide>(group, vector<double> {}));

  const int neutron = ParticleType::neutron().transport_index();
  data::energy_min[neutron] = 1.0;
  data::energy_max[neutron] = 5.0;
  data::nuclides[0]->init_grid();
}

Particle make_particle(double sqrtkT)
{
  Particle p;
  p.type() = ParticleType::neutron();
  p.E() = 2.5;
  p.E_last() = 2.5;
  p.sqrtkT() = sqrtkT;
  p.sqrtkT_last() = 0.0;
  p.stream() = STREAM_TRACKING;
  for (int i = 0; i < N_STREAMS; ++i) {
    p.seeds(i) = 17u + static_cast<uint64_t>(i);
  }
  return p;
}

NuclideMicroXS calculate_at(double sqrtkT, uint64_t* seed = nullptr)
{
  auto& nuc = *data::nuclides[0];
  auto p = make_particle(sqrtkT);
  if (seed)
    *seed = p.seeds(STREAM_TRACKING);
  nuc.calculate_xs(C_NONE, 0, 0.0, p);
  if (seed)
    *seed = p.seeds(STREAM_TRACKING);
  return p.neutron_xs(0);
}

} // namespace

TEST_CASE("nuclide temperature interpolation handles the high endpoint")
{
  GlobalState globals;
  load_temperature_fixture({0.015625, 0.0625});

  SECTION("below the available range snaps to the first table")
  {
    const auto& xs = calculate_at(std::nextafter(0.125, 0.0));
    REQUIRE(xs.index_temp == 0);
    REQUIRE(xs.absorption == Catch::Approx(10.0));
  }

  SECTION("the exact high endpoint snaps to the last table")
  {
    uint64_t seed;
    const auto& xs = calculate_at(0.25, &seed);
    REQUIRE(xs.index_temp == 1);
    REQUIRE(xs.absorption == Catch::Approx(11.0));
    REQUIRE(seed == 17u);
  }

  SECTION("above the available range snaps to the last table")
  {
    const auto& xs = calculate_at(std::nextafter(0.25, 1.0));
    REQUIRE(xs.index_temp == 1);
    REQUIRE(xs.absorption == Catch::Approx(11.0));
  }
}

TEST_CASE("nuclide temperature interpolation keeps interior sampling behavior")
{
  GlobalState globals;
  load_temperature_fixture({0.015625, 0.03515625, 0.0625});

  uint64_t seed;
  const auto& xs = calculate_at(0.1875, &seed);
  REQUIRE(xs.index_temp >= 1);
  REQUIRE(xs.index_temp <= 2);
  REQUIRE(seed != 17u);
}

TEST_CASE("single-temperature data still reverts interpolation to nearest")
{
  GlobalState globals;
  load_temperature_fixture({0.0625});

  REQUIRE(settings::temperature_method == TemperatureMethod::NEAREST);
  const auto& xs = calculate_at(0.25);
  REQUIRE(xs.index_temp == 0);
  REQUIRE(xs.absorption == Catch::Approx(10.0));
}
