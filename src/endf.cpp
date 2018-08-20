#include "openmc/endf.h"

#include <algorithm> // for copy
#include <cmath>     // for log, exp
#include <iterator>  // for back_inserter

#include "xtensor/xarray.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/search.h"

namespace openmc {

//==============================================================================
// Functions
//==============================================================================

Interpolation int2interp(int i)
{
  switch (i) {
  case 1:
    return Interpolation::histogram;
  case 2:
    return Interpolation::lin_lin;
  case 3:
    return Interpolation::lin_log;
  case 4:
    return Interpolation::log_lin;
  case 5:
    return Interpolation::log_log;
  }
}

bool is_fission(int mt)
{
  return mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38;
}

//==============================================================================
// Polynomial implementation
//==============================================================================

Polynomial::Polynomial(hid_t dset)
{
  // Read coefficients into a vector
  read_dataset(dset, coef_);
}

double Polynomial::operator()(double x) const
{
  // Use Horner's rule to evaluate polynomial. Note that coefficients are
  // ordered in increasing powers of x.
  double y = 0.0;
  for (auto c = coef_.crbegin(); c != coef_.crend(); ++c) {
    y = y*x + *c;
  }
  return y;
}

//==============================================================================
// Tabulated1D implementation
//==============================================================================

Tabulated1D::Tabulated1D(hid_t dset)
{
  read_attribute(dset, "breakpoints", nbt_);
  n_regions_ = nbt_.size();

  // Change 1-indexing to 0-indexing
  for (auto& b : nbt_) --b;

  std::vector<int> int_temp;
  read_attribute(dset, "interpolation", int_temp);

  // Convert vector of ints into Interpolation
  for (const auto i : int_temp)
    int_.push_back(int2interp(i));

  xt::xarray<double> arr;
  read_dataset(dset, arr);

  auto xs = xt::view(arr, 0);
  auto ys = xt::view(arr, 1);

  std::copy(xs.begin(), xs.end(), std::back_inserter(x_));
  std::copy(ys.begin(), ys.end(), std::back_inserter(y_));
  n_pairs_ = x_.size();
}

double Tabulated1D::operator()(double x) const
{
  // find which bin the abscissa is in -- if the abscissa is outside the
  // tabulated range, the first or last point is chosen, i.e. no interpolation
  // is done outside the energy range
  int i;
  if (x < x_[0]) {
    return y_[0];
  } else if (x > x_[n_pairs_ - 1]) {
    return y_[n_pairs_ - 1];
  } else {
    i = lower_bound_index(x_.begin(), x_.end(), x);
  }

  // determine interpolation scheme
  Interpolation interp;
  if (n_regions_ == 0) {
    interp = Interpolation::lin_lin;
  } else if (n_regions_ == 1) {
    interp = int_[0];
  } else if (n_regions_ > 1) {
    for (int j = 0; j < n_regions_; ++j) {
      if (i < nbt_[j]) {
        interp = int_[j];
        break;
      }
    }
  }

  // handle special case of histogram interpolation
  if (interp == Interpolation::histogram) return y_[i];

  // determine bounding values
  double x0 = x_[i];
  double x1 = x_[i + 1];
  double y0 = y_[i];
  double y1 = y_[i + 1];

  // determine interpolation factor and interpolated value
  double r;
  switch (interp) {
  case Interpolation::lin_lin:
    r = (x - x0)/(x1 - x0);
    return y0 + r*(y1 - y0);
  case Interpolation::lin_log:
    r = log(x/x0)/log(x1/x0);
    return y0 + r*(y1 - y0);
  case Interpolation::log_lin:
    r = (x - x0)/(x1 - x0);
    return y0*exp(r*log(y1/y0));
  case Interpolation::log_log:
    r = log(x/x0)/log(x1/x0);
    return y0*exp(r*log(y1/y0));
  }
}

//==============================================================================
// CoherentElasticXS implementation
//==============================================================================

CoherentElasticXS::CoherentElasticXS(hid_t dset)
{
  // Read 2D array from dataset
  xt::xarray<double> arr;
  read_dataset(dset, arr);

  // Get views for Bragg edges and structure factors
  auto E = xt::view(arr, 0);
  auto s = xt::view(arr, 1);

  // Copy Bragg edges and partial sums of structure factors
  std::copy(E.begin(), E.end(), std::back_inserter(bragg_edges_));
  std::copy(s.begin(), s.end(), std::back_inserter(factors_));
}

double CoherentElasticXS::operator()(double E) const
{
  if (E < bragg_edges_[0]) {
    // If energy is below that of the lowest Bragg peak, the elastic cross
    // section will be zero
    return 0.0;
  } else {
    auto i_grid = lower_bound_index(bragg_edges_.begin(), bragg_edges_.end(), E);
    return factors_[i_grid] / E;
  }
}

} // namespace openmc
