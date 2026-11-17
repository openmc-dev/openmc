#include "openmc/secondary_data.h"

#include <algorithm>
#include <cmath>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "openmc/chain.h"
#include "openmc/distribution.h"
#include "openmc/hdf5_interface.h"
#include "openmc/reaction_product.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/simulation.h"

namespace {
using namespace openmc;

struct File {
  hid_t id;
  File()
  {
    hid_t properties = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_core(properties, 4096, false);
    id = H5Fcreate("secondary-test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, properties);
    H5Pclose(properties);
  }
  ~File() { H5Fclose(id); }
};

size_t field(const SecondarySnapshot& s, size_t parent, const std::string& name)
{
  for (auto id : s.nodes[parent].children)
    if (s.nodes[id].name == name)
      return id;
  FAIL("Missing snapshot field: " + name);
  return 0;
}

void f1d(hid_t parent, const char* name, double value)
{
  auto data = tensor::zeros<double>({2, 3});
  for (int i = 0; i < 3; ++i) {
    data(0, i) = i + 1.0;
    data(1, i) = value * (i + 1.0);
  }
  write_dataset(parent, name, data);
  auto d = open_dataset(parent, name);
  write_attribute(d, "type", "Tabulated1D");
  write_attribute(d, "breakpoints", vector<int> {2, 3});
  write_attribute(d, "interpolation", vector<int> {1, 2});
  close_dataset(d);
}

// Exercise the actual HDF5 constructors, including laws absent from evaluated
// data installed on a test machine. No test-only native distribution subclass.
ReactionProduct product(
  hid_t file, const std::string& type, int angle_interp = 0)
{
  write_attribute(file, "particle", "neutron");
  write_attribute(file, "emission_mode", "prompt");
  write_attribute(file, "n_distribution", 1);
  write_dataset(file, "yield", vector<double> {1.0});
  auto y = open_dataset(file, "yield");
  write_attribute(y, "type", "Polynomial");
  close_dataset(y);
  auto d = create_group(file, "distribution_0");
  bool joint =
    type == "correlated" || type == "kalbach-mann" || type == "nbody";
  write_attribute(d, "type", joint ? type : std::string("uncorrelated"));
  auto e = joint || type == "absent" ? d : create_group(d, "energy");
  if (!joint && type != "absent")
    write_attribute(e, "type", type);
  if (type == "discrete_photon") {
    write_attribute(e, "primary_flag", 2);
    write_attribute(e, "energy", 0.5);
    write_attribute(e, "atomic_weight_ratio", 4.0);
  } else if (type == "level") {
    write_attribute(e, "threshold", 0.5);
    write_attribute(e, "mass_ratio", 0.8);
  } else if (type == "maxwell" || type == "evaporation" || type == "watt") {
    write_attribute(e, "u", -1.0);
    if (type == "watt") {
      f1d(e, "a", 0.5);
      f1d(e, "b", 0.25);
    } else {
      f1d(e, "theta", 0.5);
    }
  } else if (type == "nbody") {
    write_attribute(e, "n_particles", 3);
    write_attribute(e, "total_mass", 3.0);
    write_attribute(e, "atomic_weight_ratio", 2.0);
    write_attribute(e, "q_value", -0.5);
  } else if (type != "absent") {
    write_dataset(e, "energy", vector<double> {1.0, 3.0});
    auto grid = open_dataset(e, "energy");
    int interp[] {2, 2};
    hsize_t shape[] {2, 1};
    write_attr(grid, 2, shape, "interpolation", H5T_NATIVE_INT, interp);
    close_dataset(grid);
    auto values = tensor::zeros<double>({joint ? 5u : 3u, 4});
    for (int i = 0; i < 4; ++i) {
      values(0, i) = i % 2;
      values(1, i) = 1.0;
      values(2, i) = i % 2;
      if (joint) {
        values(3, i) = type == "correlated" ? 2.0 : 0.4;
        values(4, i) = type == "correlated" ? 2.0 * i : 1.2;
      }
    }
    const char* name = type == "correlated" ? "energy_out" : "distribution";
    write_dataset(e, name, values);
    auto table = open_dataset(e, name);
    write_attribute(table, "offsets", vector<int> {0, 2});
    write_attribute(table, "interpolation", vector<int> {2, 2});
    write_attribute(table, "n_discrete_lines", vector<int> {0, 0});
    close_dataset(table);
    if (type == "correlated") {
      auto mu = tensor::zeros<double>({3, 8});
      for (int i = 0; i < 8; ++i) {
        mu(0, i) = i % 2 ? 1.0 : -1.0;
        mu(1, i) = 0.5;
        mu(2, i) = i % 2;
      }
      write_dataset(e, "mu", mu);
    }
  }
  if (angle_interp) {
    auto a = create_group(d, "angle");
    write_dataset(a, "energy", vector<double> {1.0, 3.0});
    auto values = tensor::zeros<double>({3, 6});
    for (int i = 0; i < 2; ++i) {
      double left = i == 0 ? 0.25 : 0.75;
      for (int j = 0; j < 3; ++j) {
        values(0, 3 * i + j) = j - 1.0;
        values(1, 3 * i + j) = angle_interp == 1
                                 ? (j == 0 ? left : 1 - left)
                                 : (i == 0 ? j * 0.5 : 1 - j * 0.5);
        values(2, 3 * i + j) = j == 0 ? 0 : (j == 2 ? 1 : left);
      }
    }
    write_dataset(a, "mu", values);
    auto mu = open_dataset(a, "mu");
    write_attribute(mu, "offsets", vector<int> {0, 3});
    write_attribute(
      mu, "interpolation", vector<int> {angle_interp, angle_interp});
    close_dataset(mu);
    close_group(a);
  }
  if (e != d)
    close_group(e);
  close_group(d);
  return ReactionProduct(file);
}
} // namespace

TEST_CASE(
  "Every native reaction law can be inspected without changing sampling")
{
  auto type =
    GENERATE(std::string("continuous"), std::string("discrete_photon"),
      std::string("level"), std::string("maxwell"), std::string("evaporation"),
      std::string("watt"), std::string("correlated"),
      std::string("kalbach-mann"), std::string("nbody"));
  CAPTURE(type);
  File f;
  auto p = product(f.id, type);
  auto& law = *p.distribution_[0];
  uint64_t seed = 17;
  double before_e, before_mu;
  law.sample(2.0, before_e, before_mu, &seed);
  auto snapshot = inspect_secondary(law);
  REQUIRE(snapshot.nodes[0].type == OPENMC_SECONDARY_OBJECT);
  auto root_type = snapshot.nodes[field(snapshot, 0, "type")].text;
  auto id = root_type == "uncorrelated" ? field(snapshot, 0, "energy") : 0;
  if (type == "watt" || type == "maxwell" || type == "evaporation") {
    auto param = field(snapshot, id, type == "watt" ? "a" : "theta");
    REQUIRE(snapshot.nodes[field(snapshot, param, "breakpoints")].integers ==
            vector<int64_t> {2, 3});
    REQUIRE(snapshot.nodes[field(snapshot, param, "interpolation")].integers ==
            vector<int64_t> {1, 2});
    REQUIRE(snapshot.nodes[field(snapshot, param, "y")].reals ==
            vector<double> {0.5, 1.0, 1.5});
    if (type == "watt") {
      auto b = field(snapshot, id, "b");
      REQUIRE(snapshot.nodes[field(snapshot, b, "y")].reals ==
              vector<double> {0.25, 0.5, 0.75});
    }
  }
  if (type == "correlated" || type == "kalbach-mann" || type == "continuous") {
    auto tables = field(snapshot, id, "tables");
    REQUIRE(snapshot.nodes[tables].children.size() == 2);
    auto table = snapshot.nodes[tables].children[0];
    REQUIRE(snapshot.nodes[field(snapshot, table, "cdf")].reals ==
            vector<double> {0.0, 1.0});
    if (type == "correlated") {
      auto angles = field(snapshot, table, "angle");
      REQUIRE(snapshot.nodes[angles].children.size() == 2);
    } else if (type == "kalbach-mann") {
      REQUIRE(snapshot.nodes[field(snapshot, table, "r")].reals ==
              vector<double> {0.4, 0.4});
      REQUIRE(snapshot.nodes[field(snapshot, table, "a")].reals ==
              vector<double> {1.2, 1.2});
    }
  }
  // Destroy the copied values; native sampling must remain bitwise identical.
  for (auto& node : snapshot.nodes) {
    for (auto& x : node.reals)
      x = -123.0;
  }
  seed = 17;
  double after_e, after_mu;
  law.sample(2.0, after_e, after_mu, &seed);
  REQUIRE(after_e == before_e);
  REQUIRE(after_mu == before_mu);
}

TEST_CASE("Absent elastic energy is represented by null")
{
  File f;
  auto p = product(f.id, "absent");
  auto s = inspect_secondary(*p.distribution_[0]);
  REQUIRE(s.nodes[field(s, 0, "energy")].type == OPENMC_SECONDARY_NULL);
  auto angle = field(s, 0, "angle");
  REQUIRE(s.nodes[field(s, angle, "type")].text == "isotropic");
}

TEST_CASE("D1S borrowed spectra preserve generic native representations")
{
  const char* xml = GENERATE(
    "<source type='discrete'><parameters>1 2 0.25 0.75</parameters></source>",
    "<source type='tabular'><parameters>0 2 0.5 0.5</parameters></source>",
    "<source type='uniform'><parameters>1 2</parameters></source>",
    "<source type='maxwell'><parameters>1</parameters></source>",
    "<source type='watt'><parameters>1 2</parameters></source>",
    "<source type='normal'><parameters>2 0.1</parameters></source>",
    "<source type='powerlaw'><parameters>1 2 1</parameters></source>",
    "<source type='mixture'><pair probability='0.3'><dist type='discrete'>"
    "<parameters>1 1</parameters></dist></pair><pair probability='0.7'>"
    "<dist type='uniform'><parameters>1 2</parameters></dist></pair></source>");
  pugi::xml_document doc;
  REQUIRE(doc.load_string(xml));
  auto spectrum = distribution_from_xml(doc.child("source"));
  DecayPhotonAngleEnergy d(spectrum.get());
  auto s = inspect_secondary(d);
  REQUIRE(s.nodes[field(s, 0, "type")].text == "decay-photon");
  REQUIRE(
    s.nodes[field(s, 0, "ownership")].text == "borrowed-from-chain-nuclide");
  auto spec = field(s, 0, "spectrum");
  REQUIRE(s.nodes[field(s, spec, "type")].text != "unsupported");
  uint64_t seed = 3, control = 3;
  double e, mu;
  d.sample(2.0, e, mu, &seed);
  REQUIRE(e == spectrum->sample(&control).first);
}

TEST_CASE("Snapshot C ABI validates selection and node bounds")
{
  REQUIRE(openmc_secondary_snapshot_create(0, -1, -1, -1, false, nullptr) ==
          OPENMC_E_INVALID_ARGUMENT);
  void* s = reinterpret_cast<void*>(1);
  REQUIRE(openmc_secondary_snapshot_create(0, -1, 0, -1, false, &s) ==
          OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(s == nullptr);
  REQUIRE(openmc_secondary_snapshot_create(-1, -1, -1, -1, false, &s) ==
          OPENMC_E_OUT_OF_BOUNDS);
  bool in_batch = simulation::in_batch;
  simulation::in_batch = true;
  int code = openmc_secondary_snapshot_create(0, -1, -1, -1, false, &s);
  uint64_t seed = 1;
  double mu;
  int sample_code =
    openmc_secondary_sample_angle(0, 2, 0, 0, 1.0, &seed, 1, &mu);
  bool initialized = simulation::initialized;
  simulation::initialized = true;
  int angle_mutation = openmc_secondary_set_angle(
    0, 2, 0, 0, nullptr, 0, nullptr, nullptr, nullptr, nullptr, 0);
  uint64_t backup = 123;
  int backup_code = openmc_nuclear_data_backup(0, 2, 0, 0, 2, &backup);
  int energy_sampling =
    openmc_secondary_sample_energy(0, 24, 0, 0, 1.0, &seed, 1, &mu);
  simulation::initialized = initialized;
  simulation::in_batch = in_batch;
  REQUIRE(angle_mutation == OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(backup_code == OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(backup == 0);
  REQUIRE(energy_sampling == OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(sample_code == OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(seed == 1);
  REQUIRE(code == OPENMC_E_INVALID_ARGUMENT);
  REQUIRE(s == nullptr);
  openmc_secondary_node result;
  REQUIRE(openmc_secondary_snapshot_node(nullptr, 0, &result) ==
          OPENMC_E_INVALID_ARGUMENT);
  SecondarySnapshot empty;
  REQUIRE(openmc_secondary_snapshot_node(&empty, 0, &result) ==
          OPENMC_E_OUT_OF_BOUNDS);
  openmc_secondary_snapshot_free(nullptr);
}

TEST_CASE("Angular snapshots match native histogram and linear sampling")
{
  int interp = GENERATE(1, 2);
  double energy = GENERATE(1.0, 2.0, 3.0);
  CAPTURE(interp, energy);
  File f;
  auto p = product(f.id, "level", interp);
  auto& law = *p.distribution_[0];
  auto s = inspect_secondary(law);
  auto angle = field(s, 0, "angle");
  REQUIRE(s.nodes[field(s, angle, "incident_energy")].reals ==
          vector<double> {1.0, 3.0});
  auto tables = field(s, angle, "tables");
  auto first = s.nodes[tables].children[0];
  REQUIRE(
    s.nodes[field(s, first, "x")].reals == vector<double> {-1.0, 0.0, 1.0});
  REQUIRE(s.nodes[field(s, first, "interpolation")].integers[0] == interp);
  REQUIRE(
    s.nodes[field(s, first, "cdf")].reals == vector<double> {0.0, 0.25, 1.0});
  REQUIRE(s.nodes[field(s, first, "pdf")].reals ==
          (interp == 1 ? vector<double> {0.25, 0.75, 0.75}
                       : vector<double> {0.0, 0.5, 1.0}));
  uint64_t seed = 42;
  constexpr int n = 100000;
  vector<int> counts(8, 0);
  for (int i = 0; i < n; ++i) {
    double e, mu;
    law.sample(energy, e, mu, &seed);
    REQUIRE(mu >= -1.0);
    REQUIRE(mu <= 1.0);
    ++counts[std::min(7, static_cast<int>((mu + 1.0) * 4))];
  }
  for (int i = 0; i < 8; ++i) {
    double lo = -1.0 + i * 0.25;
    double hi = lo + 0.25;
    // Integrate the exposed first and second PDFs; incident stochastic mixing
    // gives a uniform PDF at the midpoint in both interpolation cases.
    double r = (energy - 1.0) / 2.0;
    double first_mass = interp == 1 ? (hi - lo) * (i < 4 ? 0.25 : 0.75)
                                    : (hi - lo) * (1 + (lo + hi) / 2) / 2;
    double mass = (1 - r) * first_mass + r * (0.25 - first_mass);
    double expected = n * mass;
    REQUIRE(std::abs(counts[i] - expected) < 6 * std::sqrt(expected));
  }
}

TEST_CASE(
  "Angular replacement preserves native interpolation and samples rebuilt CDF")
{
  int interp = GENERATE(1, 2);
  File f;
  auto p = product(f.id, "level", interp);
  auto& law = dynamic_cast<UncorrelatedAngleEnergy&>(*p.distribution_[0]);
  vector<double> energy {1.0, 3.0};
  vector<int> offsets {0, 3, 6}, interpolation {interp, interp};
  vector<double> x {-1, 0, 1, -1, 0, 1}, pdf {1, 2, 3, 1, 2, 3};
  replace_secondary_angle(law.angle(), energy, offsets, interpolation, x, pdf);
  auto s = inspect_secondary(law);
  auto angle = field(s, 0, "angle");
  auto tables = field(s, angle, "tables");
  auto first = s.nodes[tables].children[0];
  auto density = s.nodes[field(s, first, "pdf")].reals;
  auto cdf = s.nodes[field(s, first, "cdf")].reals;
  double norm = interp == 1 ? 3.0 : 4.0;
  REQUIRE(density == vector<double> {1 / norm, 2 / norm, 3 / norm});
  REQUIRE(std::abs(cdf[1] - (interp == 1 ? 1 / 3.0 : .375)) < 1.e-14);
  REQUIRE(std::abs(cdf[2] - 1) < 1.e-14);
  uint64_t seed = 31;
  int negative = 0;
  for (int i = 0; i < 100000; ++i) {
    double outgoing, mu;
    law.sample(2.0, outgoing, mu, &seed);
    if (mu < 0)
      ++negative;
  }
  REQUIRE(std::abs(negative / 100000.0 - cdf[1]) < .01);
  pdf[0] = 100;
  pdf[5] = -1;
  REQUIRE_THROWS_AS(replace_secondary_angle(
                      law.angle(), energy, offsets, interpolation, x, pdf),
    std::invalid_argument);
  auto after = inspect_secondary(law);
  REQUIRE(after.nodes[field(after, first, "pdf")].reals == density);
}
