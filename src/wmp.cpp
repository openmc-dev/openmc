#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"

#include <fmt/core.h>

#include <cmath>

namespace openmc {

//========================================================================
// WindowedeMultipole implementation
//========================================================================

WindowedMultipole::WindowedMultipole(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  // Read scalar values.
  read_dataset(group, "spacing", spacing_);
  read_dataset(group, "sqrtAWR", sqrt_awr_);
  read_dataset(group, "E_min", E_min_);
  read_dataset(group, "E_max", E_max_);

  // Read the "data" array.  Use its shape to figure out the number of poles
  // and residue types in this data.
  read_dataset(group, "data", data_);
  int n_residues = data_.shape()[1] - 1;

  // Check to see if this data includes fission residues.
  fissionable_ = (n_residues == 3);

  // Read the "windows" array and use its shape to figure out the number of
  // windows.
  read_dataset(group, "windows", windows_);
  int n_windows = windows_.shape()[0];

  // Read the "broaden_poly" arrays.
  read_dataset(group, "broaden_poly", broaden_poly_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }

  // Read the "curvefit" array.
  read_dataset(group, "curvefit", curvefit_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }
  fit_order_ = curvefit_.shape()[1] - 1;
}

std::tuple<double, double, double>
WindowedMultipole::evaluate(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = (sqrtE - std::sqrt(E_min_)) / spacing_;
  int startw = windows_(i_window, 0) - 1;
  int endw = windows_(i_window, 1) - 1;

  // Initialize the ouptut cross sections
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && broaden_poly_(i_window)) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::vector<double> broadened_polynomials(fit_order_ + 1);
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly, FIT_S) * broadened_polynomials[i_poly];
      sig_a += curvefit_(i_window, i_poly, FIT_A) * broadened_polynomials[i_poly];
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * broadened_polynomials[i_poly];
      }
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly, FIT_S) * temp;
      sig_a += curvefit_(i_window, i_poly, FIT_A) * temp;
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * temp;
      }
      temp *= sqrtE;
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> psi_chi = -1.0i / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi / E;
      sig_s += (data_(i_pole, MP_RS) * c_temp).real();
      sig_a += (data_(i_pole, MP_RA) * c_temp).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    if (endw >= startw) {
      for (int i_pole = startw; i_pole <= endw; ++i_pole) {
        std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
        std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
        sig_s += (data_(i_pole, MP_RS) * w_val).real();
        sig_a += (data_(i_pole, MP_RA) * w_val).real();
        if (fissionable_) {
          sig_f += (data_(i_pole, MP_RF) * w_val).real();
        }
      }
    }
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

std::tuple<double, double, double>
WindowedMultipole::evaluate_deriv(double E, double sqrtkT)
{
  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;
  double T = sqrtkT*sqrtkT / K_BOLTZMANN;

  if (sqrtkT == 0.0) {
    fatal_error("Windowed multipole temperature derivatives are not implemented"
      " for 0 Kelvin cross sections.");
  }

  // Locate us
  int i_window = (sqrtE - std::sqrt(E_min_)) / spacing_;
  int startw = windows_(i_window, 0) - 1;
  int endw = windows_(i_window, 1) - 1;

  // Initialize the ouptut cross sections.
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // TODO Polynomials: Some of the curvefit polynomials Doppler broaden so
  // rigorously we should be computing the derivative of those.  But in
  // practice, those derivatives are only large at very low energy and they
  // have no effect on reactor calculations.

  // ==========================================================================
  // Add the contribution from the poles in this window.

  double dopp = sqrt_awr_ / sqrtkT;
  if (endw >= startw) {
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
      std::complex<double> w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
      sig_s += (data_(i_pole, MP_RS) * w_val).real();
      sig_a += (data_(i_pole, MP_RA) * w_val).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * w_val).real();
      }
    }
    sig_s *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
    sig_a *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
    sig_f *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

//========================================================================
// Non-member functions
//========================================================================

void check_wmp_version(hid_t file)
{
  if (attribute_exists(file, "version")) {
    std::array<int, 2> version;
    read_attribute(file, "version", version);
    if (version[0] != WMP_VERSION[0]) {
      fatal_error(fmt::format(
        "WMP data format uses version {}.{} whereas your installation of "
        "OpenMC expects version {}.x data.",
        version[0], version[1], WMP_VERSION[0]));
    }
  } else {
    fatal_error(fmt::format("WMP data does not indicate a version. Your "
      "installation of OpenMC expects version {}x data.", WMP_VERSION[0]));
  }
}

void read_multipole_data(int i_nuclide)
{
  // Look for WMP data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::wmp, nuc->name_});

  // If no WMP library for this nuclide, just return
  if (it == data::library_map.end()) return;

  // Check if WMP library exists
  int idx = it->second;
  std::string& filename = data::libraries[idx].path_;

  // Display message
  write_message("Reading " + nuc->name_ + " WMP data from " + filename, 6);

  // Open file and make sure version is sufficient
  hid_t file = file_open(filename, 'r');
  check_wmp_version(file);

  // Read nuclide data from HDF5
  hid_t group = open_group(file, nuc->name_.c_str());
  nuc->multipole_ = std::make_unique<WindowedMultipole>(group);
  close_group(group);
  file_close(file);
}

} // namespace openmc
