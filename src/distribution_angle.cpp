#include "openmc/distribution_angle.h"

#include <cmath>  // for abs, copysign
#include <vector> // for vector

#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

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
  std::vector<int> offsets;
  std::vector<int> interp;
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
      n = offsets[i+1] - j;
    } else {
      n = temp.shape()[1] - j;
    }

    // Create and initialize tabular distribution
    auto xs = xt::view(temp, 0, xt::range(j, j+n));
    auto ps = xt::view(temp, 1, xt::range(j, j+n));
    auto cs = xt::view(temp, 2, xt::range(j, j+n));
    std::vector<double> x {xs.begin(), xs.end()};
    std::vector<double> p {ps.begin(), ps.end()};
    std::vector<double> c {cs.begin(), cs.end()};

    // To get answers that match ACE data, for now we still use the tabulated
    // CDF values that were passed through to the HDF5 library. At a later
    // time, we can remove the CDF values from the HDF5 library and
    // reconstruct them using the PDF
    Tabular* mudist = new Tabular{x.data(), p.data(), n, int2interp(interp[i]),
                                  c.data()};

    distribution_.emplace_back(mudist);
  }
}

double AngleDistribution::sample(double E, uint64_t* seed) const
{
  // Determine number of incoming energies
  auto n = energy_.size();

  // Find energy bin and calculate interpolation factor -- if the energy is
  // outside the range of the tabulated energies, choose the first or last bins
  int i;
  double r;
  if (E < energy_[0]) {
    i = 0;
    r = 0.0;
  } else if (E > energy_[n - 1]) {
    i = n - 2;
    r = 1.0;
  } else {
    i = lower_bound_index(energy_.begin(), energy_.end(), E);
    r = (E - energy_[i])/(energy_[i+1] - energy_[i]);
  }

  // Sample between the ith and (i+1)th bin
  if (r > prn(seed)) ++i;

  // Sample i-th distribution
  double mu = distribution_[i]->sample(seed);

  // Make sure mu is in range [-1,1] and return
  if (std::abs(mu) > 1.0) mu = std::copysign(1.0, mu);
  return mu;
}

} // namespace openmc
