#include "openmc/distribution_angle.h"

#include <cmath> // for abs, copysign

#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/vector.h" // for vector

namespace openmc {

//==============================================================================
// AngleDistribution implementation
//==============================================================================

AngleDistribution::AngleDistribution(hid_t group)
{
  // Get incoming energies
  read_dataset(group, "energy", energy_);
  int n_energy = energy_.size();

  // Get outgoing energy distribution data
  vector<int> offsets;
  vector<int> interp;
  hid_t dset = open_dataset(group, "mu");
  read_attribute(dset, "offsets", offsets);
  read_attribute(dset, "interpolation", interp);
  xt::xarray<double> temp;
  read_dataset(dset, temp);
  close_dataset(dset);

  for (int i = 0; i < n_energy; ++i) {
    // Determine number of outgoing energies
    int j = offsets[i];
    int n;
    if (i < n_energy - 1) {
      n = offsets[i + 1] - j;
    } else {
      n = temp.shape()[1] - j;
    }

    // Create and initialize tabular distribution
    auto xs = xt::view(temp, 0, xt::range(j, j + n));
    auto ps = xt::view(temp, 1, xt::range(j, j + n));
    auto cs = xt::view(temp, 2, xt::range(j, j + n));
    vector<double> x {xs.begin(), xs.end()};
    vector<double> p {ps.begin(), ps.end()};
    vector<double> c {cs.begin(), cs.end()};

    // To get answers that match ACE data, for now we still use the tabulated
    // CDF values that were passed through to the HDF5 library. At a later
    // time, we can remove the CDF values from the HDF5 library and
    // reconstruct them using the PDF
    Tabular* mudist =
      new Tabular {x.data(), p.data(), n, int2interp(interp[i]), c.data()};

    distribution_.emplace_back(mudist);
  }
}

double AngleDistribution::sample(double E, uint64_t* seed) const
{
  // Find energy bin and calculate interpolation factor
  int i;
  double r;
  get_energy_index(energy_, E, i, r);

  // Sample between the ith and (i+1)th bin
  if (r > prn(seed))
    ++i;

  // Sample i-th distribution
  double mu = distribution_[i]->sample(seed);

  // Make sure mu is in range [-1,1] and return
  if (std::abs(mu) > 1.0)
    mu = std::copysign(1.0, mu);
  return mu;
}

double AngleDistribution::evaluate(double E, double mu) const
{
  // Find energy bin and calculate interpolation factor
  int i;
  double r;
  get_energy_index(energy_, E, i, r);

  double pdf = 0.0;
  if (r > 0.0)
    pdf += r * distribution_[i + 1]->evaluate(mu);
  if (r < 1.0)
    pdf += (1.0 - r) * distribution_[i]->evaluate(mu);
  return pdf;
}

} // namespace openmc
