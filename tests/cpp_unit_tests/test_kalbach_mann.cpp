#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/hdf5_interface.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/tensor.h"
#include "openmc/vector.h"

using namespace openmc;
using Catch::Matchers::WithinAbs;

namespace {

// A Kalbach-Mann distribution with constant parameters R and A at two incident
// energies, and outgoing energies uniform on [0, 1 MeV]
KalbachMann make_kalbach_mann(
  const std::string& filename, const char* particle, double R, double A)
{
  hid_t file = file_open(filename, 'w');
  hid_t group = create_group(file, "km");
  if (particle)
    write_attribute(group, "particle", particle);

  write_dataset(group, "energy", vector<double> {1.0e7, 2.0e7});
  hid_t dset = open_dataset(group, "energy");
  // One region with histogram interpolation: breakpoints, then laws
  vector<int> interpolation {2, 1};
  hsize_t dims[] {2, 1};
  write_attr(
    dset, 2, dims, "interpolation", H5T_NATIVE_INT, interpolation.data());
  close_dataset(dset);

  // Rows are outgoing energy, PDF, CDF, R and A
  auto eout = tensor::zeros<double>({5, 4});
  for (int j = 0; j < 4; ++j) {
    bool upper = (j % 2 == 1);
    eout(0, j) = upper ? 1.0e6 : 0.0;
    eout(1, j) = upper ? 0.0 : 1.0e-6;
    eout(2, j) = upper ? 1.0 : 0.0;
    eout(3, j) = R;
    eout(4, j) = A;
  }
  write_dataset(group, "distribution", eout);
  dset = open_dataset(group, "distribution");
  write_attribute(dset, "offsets", vector<int> {0, 2});
  write_attribute(dset, "interpolation", vector<int> {1, 1});
  write_attribute(dset, "n_discrete_lines", vector<int> {0, 0});
  close_dataset(dset);

  KalbachMann dist {group};
  close_group(group);
  file_close(file);
  std::remove(filename.c_str());
  return dist;
}

// Check that the angular PDF is normalized and agrees with the cosines that
// sample() produces
void check_kalbach_mann_pdf(const KalbachMann& dist)
{
  const double E_in = 1.5e7;
  uint64_t seed = 1;
  double E_out;

  // Integrate the PDF over each of n_bins equal cosine bins with Simpson's rule
  const int n_bins = 20;
  const int n_sub = 20;
  const double width = 2.0 / n_bins;
  vector<double> prob(n_bins, 0.0);
  for (int b = 0; b < n_bins; ++b) {
    double h = width / n_sub;
    for (int k = 0; k <= n_sub; ++k) {
      double mu = -1.0 + b * width + k * h;
      double w = (k == 0 || k == n_sub) ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0);
      prob[b] +=
        w * h / 3.0 * dist.sample_energy_and_pdf(E_in, mu, E_out, &seed);
    }
  }
  double total = 0.0;
  for (double p : prob)
    total += p;
  REQUIRE_THAT(total, WithinAbs(1.0, 1e-8));

  const int n = 400000;
  vector<int> counts(n_bins, 0);
  for (int i = 0; i < n; ++i) {
    double mu;
    dist.sample(E_in, E_out, mu, &seed);
    int b = std::min(static_cast<int>((mu + 1.0) / width), n_bins - 1);
    ++counts[b];
  }
  for (int b = 0; b < n_bins; ++b) {
    double expected = n * prob[b];
    double sigma = std::sqrt(expected * (1.0 - prob[b]));
    INFO("bin " << b << ": expected " << expected << ", sampled " << counts[b]);
    REQUIRE(std::abs(counts[b] - expected) < 5.0 * sigma);
  }
}

} // namespace

// For an incident photon the Kalbach-Mann angular distribution is isotropic
// with probability 1 - R and proportional to exp(A mu) with probability R,
// rather than the symmetric form used for an incident neutron. The PDF has to
// describe the distribution that sample() draws from.
TEST_CASE("Kalbach-Mann angular PDF matches sampling")
{
  SECTION("Incident neutron")
  {
    check_kalbach_mann_pdf(
      make_kalbach_mann("test_km_neutron.h5", nullptr, 0.4, 1.5));
  }
  SECTION("Incident photon")
  {
    check_kalbach_mann_pdf(
      make_kalbach_mann("test_km_photon.h5", "photon", 0.4, 1.5));
  }
}
