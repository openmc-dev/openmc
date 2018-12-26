#ifndef OPENMC_WMP_H
#define OPENMC_WMP_H

#include "hdf5.h"
#include "xtensor/xtensor.hpp"

#include <string>

namespace openmc {

class WindowedMultipole {
public:
  // Types, aliases
  using cdouble = std::complex<double>;

  // Constructors, destructors
  WindowedMultipole(hid_t group);

  // Data members
  std::string name_; //!< Name of nuclide
  bool fissionable_; //!< Is the nuclide fissionable?
  xt::xtensor<cdouble, 2> data_; //!< Poles and residues
  double sqrt_awr_; //!< Square root of atomic weight ratio
  double E_min_; //!< Minimum energy in [eV]
  double E_max_; //!< Maximum energy in [eV]
  double spacing_; //!< Spacing in sqrt(E) space
  int fit_order_; //!< Order of the fit
  xt::xtensor<int, 2> windows_; //!< Indices of pole at start/end of window
  xt::xtensor<double, 2> curvefit_; //!< Fitting function (reaction, coeff index, window index)
  xt::xtensor<bool, 1> broaden_poly_; //!< Whether to broaden curvefit
};

} // namespace openmc

#endif // OPENMC_WMP_H
