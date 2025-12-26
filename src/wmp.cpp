#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/error.h" // for writing messages
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"

#include <fmt/core.h>

#include <algorithm> // for min
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
  read_dataset(group, "spacing", inv_spacing_);
  inv_spacing_ = 1.0 / inv_spacing_;
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
  xt::xtensor<int, 2> windows;
  read_dataset(group, "windows", windows);
  int n_windows = windows.shape()[0];
  windows -= 1; // Adjust to 0-based indices

  // Read the "broaden_poly" arrays.
  xt::xtensor<bool, 1> broaden_poly;
  read_dataset(group, "broaden_poly", broaden_poly);
  if (n_windows != broaden_poly.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
                "array shape in WMP library for " +
                name_ + ".");
  }

  // Read the "curvefit" array.
  read_dataset(group, "curvefit", curvefit_);
  if (n_windows != curvefit_.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
                "array shape in WMP library for " +
                name_ + ".");
  }
  fit_order_ = curvefit_.shape()[1] - 1;

  // Check the code is compiling to work with sufficiently high fit order
  if (fit_order_ + 1 > MAX_POLY_COEFFICIENTS) {
    fatal_error(fmt::format(
      "Need to compile with WindowedMultipole::MAX_POLY_COEFFICIENTS = {}",
      fit_order_ + 1));
  }

  // Move window information into a vector
  window_info_.resize(n_windows);
  for (int i = 0; i < n_windows; ++i) {
    window_info_[i].index_start = windows(i, 0);
    window_info_[i].index_end = windows(i, 1);
    window_info_[i].broaden_poly = broaden_poly[i];
  }
}

std::tuple<double, double, double> WindowedMultipole::evaluate(
  double E, double sqrtkT) const
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Initialize the ouptut cross sections
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    array<double, MAX_POLY_COEFFICIENTS> broadened_polynomials;
    broaden_wmp_polynomials(
      E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s +=
        curvefit_(i_window, i_poly, FIT_S) * broadened_polynomials[i_poly];
      sig_a +=
        curvefit_(i_window, i_poly, FIT_A) * broadened_polynomials[i_poly];
      if (fissionable_) {
        sig_f +=
          curvefit_(i_window, i_poly, FIT_F) * broadened_polynomials[i_poly];
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
    for (int i_pole = window.index_start; i_pole <= window.index_end;
         ++i_pole) {
      std::complex<double> psi_chi = -1.0i / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi * invE;
      sig_s += (data_(i_pole, MP_RS) * c_temp).real();
      sig_a += (data_(i_pole, MP_RA) * c_temp).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    for (int i_pole = window.index_start; i_pole <= window.index_end;
         ++i_pole) {
      std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
      std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
      sig_s += (data_(i_pole, MP_RS) * w_val).real();
      sig_a += (data_(i_pole, MP_RA) * w_val).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * w_val).real();
      }
    }
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

std::tuple<double, double, double> WindowedMultipole::evaluate_deriv(
  double E, double sqrtkT) const
{
  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;
  double T = sqrtkT * sqrtkT / K_BOLTZMANN;

  if (sqrtkT == 0.0) {
    fatal_error("Windowed multipole temperature derivatives are not implemented"
                " for 0 Kelvin cross sections.");
  }

  // Locate us
  int i_window = (sqrtE - std::sqrt(E_min_)) * inv_spacing_;
  const auto& window {window_info_[i_window]};

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
  for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
    std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
    std::complex<double> w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
    sig_s += (data_(i_pole, MP_RS) * w_val).real();
    sig_a += (data_(i_pole, MP_RA) * w_val).real();
    if (fissionable_) {
      sig_f += (data_(i_pole, MP_RF) * w_val).real();
    }
  }
  double norm = -0.5 * sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  sig_s *= norm;
  sig_a *= norm;
  sig_f *= norm;

  return std::make_tuple(sig_s, sig_a, sig_f);
}

//========================================================================
// WindowedMultipole parameter derivatives
//========================================================================

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_pole_deriv_total(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = 2*(window.index_end - window.index_start + 1)*data_.shape()[1];
  std::vector<double> derivative(size, 0.0);
  int vector_start = window.index_start * 2 * data_.shape()[1];

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
      int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
      std::complex<double> psi_chi = 1.0 / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi / E;

      derivative[start_idx]   = (-1.0 * (data_(i_pole, MP_RS) + data_(i_pole, MP_RA)) * psi_chi * c_temp).imag();
      derivative[start_idx+1] = (-1.0 * (data_(i_pole, MP_RS) + data_(i_pole, MP_RA)) * psi_chi * c_temp).real();
      derivative[start_idx+2] = (c_temp).imag();
      derivative[start_idx+3] = (c_temp).real();
      derivative[start_idx+4] = (c_temp).imag();
      derivative[start_idx+5] = (c_temp).real();
      // Total = Scatter + Absorption, therefore changing r_f does not affect Total Cross Section
      //if (fissionable_){
      //  derivative[start_idx]   += (-1.0 * data_(i_pole, MP_RF) * psi_chi * c_temp).imag();
      //  derivative[start_idx+1] += (-1.0 * data_(i_pole, MP_RF) * psi_chi * c_temp).real();
      //  derivative[start_idx+6] = (c_temp).imag();
      //  derivative[start_idx+7] = (c_temp).real();
      //}
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    if (window.index_end >= window.index_start) {
      for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
        int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
        std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
        std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
        std::complex<double> dw_val = w_derivative(z,1) * dopp * invE * SQRT_PI;

        derivative[start_idx]   = (-1.0 * dopp * (data_(i_pole, MP_RS) + data_(i_pole, MP_RA)) * dw_val).real();
        derivative[start_idx+1] = (-1.0i * dopp * (data_(i_pole, MP_RS) + data_(i_pole, MP_RA)) * dw_val).real();
        derivative[start_idx+2] = (w_val).real();
        derivative[start_idx+3] = (1.0i * w_val).real();
        derivative[start_idx+4] = (w_val).real();
        derivative[start_idx+5] = (1.0i * w_val).real();
        // Total = Scatter + Absorption, therefore changing r_f does not affect Total Cross Section
        //if (fissionable_){
        //  derivative[start_idx]   += (-1.0 * data_(i_pole, MP_RF) * dopp * dw_val).real();
        //  derivative[start_idx+1] += (-1.0i * data_(i_pole, MP_RF) * dopp * dw_val).real();
        //  derivative[start_idx+6] = (w_val).real();
        //  derivative[start_idx+7] = (1.0i * w_val).real();
        //}
      }
    }
  }

  return std::make_pair(vector_start, derivative);
}

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_pole_deriv_scatter(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = 2*(window.index_end - window.index_start + 1)*data_.shape()[1];
  std::vector<double> derivative(size, 0.0);
  int vector_start = window.index_start * 2 * data_.shape()[1];

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
      int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
      std::complex<double> psi_chi = 1.0 / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi / E;

      derivative[start_idx]   = (-1.0 * data_(i_pole, MP_RS) * psi_chi * c_temp).imag();
      derivative[start_idx+1] = (-1.0 * data_(i_pole, MP_RS) * psi_chi * c_temp).real();
      derivative[start_idx+2] = (c_temp).imag();
      derivative[start_idx+3] = (c_temp).real();
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    if (window.index_end >= window.index_start) {
      for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
        int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
        std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
        std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
        std::complex<double> dw_val = w_derivative(z,1) * dopp * invE * SQRT_PI;

        derivative[start_idx]   = (-1.0 * dopp * data_(i_pole, MP_RS) * dw_val).real();
        derivative[start_idx+1] = (-1.0i * dopp * data_(i_pole, MP_RS) * dw_val).real();
        derivative[start_idx+2] = (w_val).real();
        derivative[start_idx+3] = (1.0i * w_val).real();
      }
    }
  }

  return std::make_pair(vector_start, derivative);
}

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_pole_deriv_fission(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = 2*(window.index_end - window.index_start + 1)*data_.shape()[1];
  std::vector<double> derivative(size, 0.0);
  int vector_start = window.index_start * 2 * data_.shape()[1];

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
      int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
      std::complex<double> psi_chi = 1.0 / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi / E;
      
      derivative[start_idx]   = (-1.0 * data_(i_pole, MP_RF) * psi_chi * c_temp).imag();
      derivative[start_idx+1] = (-1.0 * data_(i_pole, MP_RF) * psi_chi * c_temp).real();
      derivative[start_idx+6] = (c_temp).imag();
      derivative[start_idx+7] = (c_temp).real();
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    if (window.index_end >= window.index_start) {
      for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
        int start_idx = (i_pole - window.index_start) * 2 * data_.shape()[1];
        std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
        std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
        std::complex<double> dw_val = w_derivative(z,1) * dopp * invE * SQRT_PI;
        
        derivative[start_idx]   = (-1.0 * data_(i_pole, MP_RF) * dopp * dw_val).real();
        derivative[start_idx+1] = (-1.0i * data_(i_pole, MP_RF) * dopp * dw_val).real();
        derivative[start_idx+6] = (w_val).real();
        derivative[start_idx+7] = (1.0i * w_val).real();
      }
    }
  }

  return std::make_pair(vector_start, derivative);
}

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_fit_deriv_total(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy 
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = curvefit_.shape()[1]*curvefit_.shape()[2];
  std::vector<double> derivative(size, 0.0);
  int vector_start = i_window*curvefit_.shape()[1]*curvefit_.shape()[2];

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::vector<double> broadened_polynomials(fit_order_ + 1);
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = broadened_polynomials[i_poly];
      derivative[i_poly + FIT_A*curvefit_.shape()[1]] = broadened_polynomials[i_poly];
      // Total = Scatter + Absorption, therefore changing r_f does not affect Total Cross Section
      //if (fissionable_) {
      //  derivative[i_poly + FIT_F*curvefit_.shape()[1]] = broadened_polynomials[i_poly];
      //}
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = temp;
      derivative[i_poly + FIT_A*curvefit_.shape()[1]] = temp;
      // Total = Scatter + Absorption, therefore changing r_f does not affect Total Cross Section
      //if (fissionable_) {
      //  derivative[i_poly + FIT_F*curvefit_.shape()[1]] = temp;
      //}
      temp *= sqrtE;
    }
  }
  
  return std::make_pair(vector_start, derivative);
}

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_fit_deriv_scatter(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = curvefit_.shape()[1];
  std::vector<double> derivative(size, 0.0);
  int vector_start = i_window*curvefit_.shape()[1]*curvefit_.shape()[2] + FIT_S*curvefit_.shape()[1];

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::vector<double> broadened_polynomials(fit_order_ + 1);
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = broadened_polynomials[i_poly];
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = temp;
      temp *= sqrtE;
    }
  }
  
  return std::make_pair(vector_start, derivative);
}

std::pair<int, std::vector<double>>
WindowedMultipole::evaluate_fit_deriv_fission(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  // int i_window = (sqrtE - std::sqrt(E_min_)) * inv_spacing_;
  
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Calculate size of vector and initialize
  int size = curvefit_.shape()[1];
  std::vector<double> derivative(size, 0.0);
  int vector_start = i_window*curvefit_.shape()[1]*curvefit_.shape()[2] + FIT_F*curvefit_.shape()[1];

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::vector<double> broadened_polynomials(fit_order_ + 1);
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = broadened_polynomials[i_poly];
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      derivative[i_poly] = temp;
      temp *= sqrtE;
    }
  }

  return std::make_pair(vector_start, derivative);
}

//========================================================================
// Non-member functions
//========================================================================

void check_wmp_version(hid_t file)
{
  if (attribute_exists(file, "version")) {
    array<int, 2> version;
    read_attribute(file, "version", version);
    if (version[0] != WMP_VERSION[0]) {
      fatal_error(fmt::format(
        "WMP data format uses version {}.{} whereas your installation of "
        "OpenMC expects version {}.x data.",
        version[0], version[1], WMP_VERSION[0]));
    }
  } else {
    fatal_error(fmt::format("WMP data does not indicate a version. Your "
                            "installation of OpenMC expects version {}x data.",
      WMP_VERSION[0]));
  }
}

void read_multipole_data(int i_nuclide)
{
  // Look for WMP data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::wmp, nuc->name_});

  // If no WMP library for this nuclide, just return
  if (it == data::library_map.end())
    return;

  // Check if WMP library exists
  int idx = it->second;
  std::string& filename = data::libraries[idx].path_;

  // Display message
  write_message(6, "Reading {} WMP data from {}", nuc->name_, filename);

  // Open file and make sure version is sufficient
  hid_t file = file_open(filename, 'r');
  check_wmp_version(file);

  // Read nuclide data from HDF5
  hid_t group = open_group(file, nuc->name_.c_str());
  nuc->multipole_ = make_unique<WindowedMultipole>(group);
  close_group(group);
  file_close(file);
}

void broaden_wmp_polynomials(double E, double dopp, int n, double factors[])
{
  // Broadening of polynomials follows procedure outlined in C. Josey, P. Ducru,
  // B. Forget, and K. Smith, "Windowed multipole for cross section Doppler
  // broadening," J. Comput. Phys., 307, 715-727 (2016).
  // https://doi.org/10.1016/j.jcp.2015.08.013

  // Factors is already pre-allocated
  double sqrtE = std::sqrt(E);
  double beta = sqrtE * dopp;
  double half_inv_dopp2 = 0.5 / (dopp * dopp);
  double quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;

  double erf_beta;    // error function of beta
  double exp_m_beta2; // exp(-beta**2)
  if (beta > 6.0) {
    // Save time, ERF(6) is 1 to machine precision.
    // beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
    erf_beta = 1.;
    exp_m_beta2 = 0.;
  } else {
    erf_beta = std::erf(beta);
    exp_m_beta2 = std::exp(-beta * beta);
  }

  // Assume that, for sure, we'll use a second order (1/E, 1/V, const)
  // fit, and no less.
  factors[0] = erf_beta / E;
  factors[1] = 1. / sqrtE;
  factors[2] =
    factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 / (beta * SQRT_PI);
  if (n > 3)
    factors[3] = factors[1] * (E + 3.0 * half_inv_dopp2);

  // Perform recursive broadening of high order components (Eq. 16)
  for (int i = 1; i < n - 3; i++) {
    double ip1_dbl = i + 1;
    factors[i + 3] =
      -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl * quarter_inv_dopp4 +
      factors[i + 1] * (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
  }
}

} // namespace openmc
