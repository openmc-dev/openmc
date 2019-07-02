#include "openmc/endf.h"

#include <algorithm> // for copy
#include <array>
#include <cmath>     // for log, exp
#include <iterator>  // for back_inserter
#include <stdexcept> // for runtime_error

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
  // TODO: We are ignoring specification of two-dimensional interpolation
  // schemes (method of corresponding points and unit base interpolation). Those
  // should be accounted for in the distribution classes somehow.

  switch (i) {
  case 1: case 11: case 21:
    return Interpolation::histogram;
  case 2: case 12: case 22:
    return Interpolation::lin_lin;
  case 3: case 13: case 23:
    return Interpolation::lin_log;
  case 4: case 14: case 24:
    return Interpolation::log_lin;
  case 5: case 15: case 25:
    return Interpolation::log_log;
  default:
    throw std::runtime_error{"Invalid interpolation code."};
  }
}

bool is_fission(int mt)
{
  return mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38;
}

bool is_disappearance(int mt)
{
  if (mt >= N_DISAPPEAR && mt <= N_DA) {
    return true;
  } else if (mt >= N_P0 && mt <= N_AC) {
    return true;
  } else if (mt == N_TA || mt == N_DT || mt == N_P3HE || mt == N_D3HE
    || mt == N_3HEA || mt == N_3P) {
    return true;
  } else {
    return false;
  }
}

bool is_inelastic_scatter(int mt)
{
  if (mt < 100) {
    if (is_fission(mt)) {
      return false;
    } else {
      return mt >= MISC && mt != 27;
    }
  } else if (mt <= 200) {
    return !is_disappearance(mt);
  } else if (mt >= N_2N0 && mt <= N_2NC) {
    return true;
  } else {
    return false;
  }
}

std::unique_ptr<Function1D>
read_function(hid_t group, const char* name)
{
  hid_t dset = open_dataset(group, name);
  std::string func_type;
  read_attribute(dset, "type", func_type);
  std::unique_ptr<Function1D> func;
  if (func_type == "Tabulated1D") {
    func = std::make_unique<Tabulated1D>(dset);
  } else if (func_type == "Polynomial") {
    func = std::make_unique<Polynomial>(dset);
  } else if (func_type == "CoherentElastic") {
    func = std::make_unique<CoherentElasticXS>(dset);
  } else if (func_type == "IncoherentElastic") {
    func = std::make_unique<IncoherentElasticXS>(dset);
  } else {
    throw std::runtime_error{"Unknown function type " + func_type +
      " for dataset " + object_name(dset)};
  }
  close_dataset(dset);
  return func;
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
  } else {
    interp = int_[0];
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
  default:
    throw std::runtime_error{"Invalid interpolation scheme."};
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

//==============================================================================
// IncoherentElasticXS implementation
//==============================================================================

IncoherentElasticXS::IncoherentElasticXS(hid_t dset)
{
  std::array<double, 2> tmp;
  read_dataset(dset, nullptr, tmp);
  bound_xs_ = tmp[0];
  debye_waller_ = tmp[1];
}

double IncoherentElasticXS::operator()(double E) const
{
  // Determine cross section using ENDF-102, Eq. (7.5)
  double W = debye_waller_;
  return bound_xs_ / 2.0 * ((1 - std::exp(-4.0*E*W))/(2.0*E*W));
}

} // namespace openmc
