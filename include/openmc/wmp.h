#ifndef OPENMC_WMP_H
#define OPENMC_WMP_H

#include "hdf5.h"
#include "xtensor/xtensor.hpp"

#include <complex>
#include <string>
#include <tuple>

#include "openmc/array.h"
#include "openmc/vector.h"

namespace openmc {

//========================================================================
// Constants
//========================================================================

// Constants that determine which value to access
constexpr int MP_EA {0}; // Pole
constexpr int MP_RS {1}; // Residue scattering
constexpr int MP_RA {2}; // Residue absorption
constexpr int MP_RF {3}; // Residue fission

// Polynomial fit indices
constexpr int FIT_S {0}; // Scattering
constexpr int FIT_A {1}; // Absorption
constexpr int FIT_F {2}; // Fission

// Multipole HDF5 file version
constexpr array<int, 2> WMP_VERSION {1, 1};

//========================================================================
// Windowed multipole data
//========================================================================

class WindowedMultipole {
public:
  // Types
  struct WindowInfo {
    int index_start;   // Index of starting pole
    int index_end;     // Index of ending pole
    bool broaden_poly; // Whether to broaden polynomial curvefit
  };

  // Constructors, destructors
  WindowedMultipole(hid_t group);

  // Methods

  //! \brief Evaluate the windowed multipole equations for cross sections in the
  //! resolved resonance regions
  //!
  //! \param E Incident neutron energy in [eV]
  //! \param sqrtkT Square root of temperature times Boltzmann constant
  //! \return Tuple of elastic scattering, absorption, and fission cross
  //! sections in [b]
  std::tuple<double, double, double> evaluate(double E, double sqrtkT) const;

  //! \brief Evaluates the windowed multipole equations for the derivative of
  //! cross sections in the resolved resonance regions with respect to
  //! temperature.
  //!
  //! \param E Incident neutron energy in [eV]
  //! \param sqrtkT Square root of temperature times Boltzmann constant
  //! \return Tuple of derivatives of elastic scattering, absorption, and
  //!         fission cross sections in [b/K]
  std::tuple<double, double, double> evaluate_deriv(
    double E, double sqrtkT) const;
    
  //! function that returns derivative std::pair<int, std::vector<double>>
  //! separate functions for poles and coefficients?
  //! separate functions for curvefit coefficients?
  //! this would make a total of 6!
  std::pair<int, std::vector<double>> evaluate_pole_deriv_total(double E, double sqrtkT);

  std::pair<int, std::vector<double>> evaluate_pole_deriv_scatter(double E, double sqrtkT);

  std::pair<int, std::vector<double>> evaluate_pole_deriv_fission(double E, double sqrtkT);

  std::pair<int, std::vector<double>> evaluate_fit_deriv_total(double E, double sqrtkT);

  std::pair<int, std::vector<double>> evaluate_fit_deriv_scatter(double E, double sqrtkT);

  std::pair<int, std::vector<double>> evaluate_fit_deriv_fission(double E, double sqrtkT);

  // Data members
  std::string name_;               //!< Name of nuclide
  double E_min_;                   //!< Minimum energy in [eV]
  double E_max_;                   //!< Maximum energy in [eV]
  double sqrt_awr_;                //!< Square root of atomic weight ratio
  double inv_spacing_;             //!< 1 / spacing in sqrt(E) space
  int fit_order_;                  //!< Order of the fit
  bool fissionable_;               //!< Is the nuclide fissionable?
  vector<WindowInfo> window_info_; // Information about a window
  xt::xtensor<double, 3>
    curvefit_; // Curve fit coefficients (window, poly order, reaction)
  xt::xtensor<std::complex<double>, 2> data_; //!< Poles and residues

  // Constant data
  static constexpr int MAX_POLY_COEFFICIENTS =
    11; //!< Max order of polynomial fit plus one
};

//========================================================================
// Non-member functions
//========================================================================

//! Check to make sure WMP library data version matches
//!
//! \param[in] file  HDF5 file object
void check_wmp_version(hid_t file);

//! \brief Checks for the existence of a multipole library in the directory and
//! loads it
//!
//! \param[in] i_nuclide  Index in global nuclides array
void read_multipole_data(int i_nuclide);

//==============================================================================
//! Doppler broadens the windowed multipole curvefit.
//!
//! The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E)...
//!
//! \param E       The energy to evaluate the broadening at
//! \param dopp    sqrt(atomic weight ratio / kT) with kT given in eV
//! \param n       The number of components to the polynomial
//! \param factors The output leading coefficient
//==============================================================================

extern "C" void broaden_wmp_polynomials(
  double E, double dopp, int n, double factors[]);

} // namespace openmc

#endif // OPENMC_WMP_H
