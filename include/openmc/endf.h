//! \file endf.h
//! Classes and functions related to the ENDF-6 format

#ifndef OPENMC_ENDF_H
#define OPENMC_ENDF_H

#include "hdf5.h"

#include "openmc/constants.h"
#include "openmc/memory.h"
#include "openmc/vector.h"

namespace openmc {

//! Convert integer representing interpolation law to enum
//! \param[in] i Intereger (e.g. 1=histogram, 2=lin-lin)
//! \return Corresponding enum value
Interpolation int2interp(int i);

//! Determine whether MT number corresponds to a fission reaction
//! \param[in] MT ENDF MT value
//! \return Whether corresponding reaction is a fission reaction
bool is_fission(int MT);

//! Determine if a given MT number is that of a disappearance reaction, i.e., a
//! reaction with no neutron in the exit channel
//! \param[in] MT ENDF MT value
//! \return Whether corresponding reaction is a disappearance reaction
bool is_disappearance(int MT);

//! Determine if a given MT number is that of an inelastic scattering reaction
//! \param[in] MT ENDF MT value
//! \return Whether corresponding reaction is an inelastic scattering reaction
bool is_inelastic_scatter(int MT);

//==============================================================================
//! Abstract one-dimensional function
//==============================================================================

class Function1D {
public:
  virtual double operator()(double x) const = 0;
  virtual ~Function1D() = default;
};

//==============================================================================
//! One-dimensional function expressed as a polynomial
//==============================================================================

class Polynomial : public Function1D {
public:
  //! Construct polynomial from HDF5 data
  //! \param[in] dset Dataset containing coefficients
  explicit Polynomial(hid_t dset);

  //! Construct polynomial from coefficients
  //! \param[in] coef Polynomial coefficients
  explicit Polynomial(vector<double> coef) : coef_(coef) {}

  //! Evaluate the polynomials
  //! \param[in] x independent variable
  //! \return Polynomial evaluated at x
  double operator()(double x) const override;

private:
  vector<double> coef_; //!< Polynomial coefficients
};

//==============================================================================
//! One-dimensional interpolable function
//==============================================================================

class Tabulated1D : public Function1D {
public:
  Tabulated1D() = default;

  //! Construct function from HDF5 data
  //! \param[in] dset Dataset containing tabulated data
  explicit Tabulated1D(hid_t dset);

  explicit Tabulated1D(vector<double>& x, vector<double>& y);

  //! Evaluate the tabulated function
  //! \param[in] x independent variable
  //! \return Function evaluated at x
  double operator()(double x) const override;

  // Accessors
  const vector<double>& x() const { return x_; }
  const vector<double>& y() const { return y_; }

private:
  std::size_t n_regions_ {0}; //!< number of interpolation regions
  vector<int> nbt_;           //!< values separating interpolation regions
  vector<Interpolation> int_; //!< interpolation schemes
  std::size_t n_pairs_;       //!< number of (x,y) pairs
  vector<double> x_;          //!< values of abscissa
  vector<double> y_;          //!< values of ordinate
};

//==============================================================================
//! Coherent elastic scattering data from a crystalline material
//==============================================================================

class CoherentElasticXS : public Function1D {
public:
  explicit CoherentElasticXS(hid_t dset);

  double operator()(double E) const override;

  const vector<double>& bragg_edges() const { return bragg_edges_; }
  const vector<double>& factors() const { return factors_; }

private:
  vector<double> bragg_edges_; //!< Bragg edges in [eV]
  vector<double> factors_;     //!< Partial sums of structure factors [eV-b]
};

//==============================================================================
//! Incoherent elastic scattering cross section
//==============================================================================

class IncoherentElasticXS : public Function1D {
public:
  explicit IncoherentElasticXS(hid_t dset);

  double operator()(double E) const override;

private:
  double bound_xs_; //!< Characteristic bound xs in [b]
  double
    debye_waller_; //!< Debye-Waller integral divided by atomic mass in [eV^-1]
};

//==============================================================================
//! Sum of multiple 1D functions
//==============================================================================

class Sum1D : public Function1D {
public:
  // Constructors
  explicit Sum1D(hid_t group);

  //! Evaluate each function and sum results
  //! \param[in] x independent variable
  //! \return Function evaluated at x
  double operator()(double E) const override;

  const unique_ptr<Function1D>& functions(int i) const { return functions_[i]; }

private:
  vector<unique_ptr<Function1D>> functions_; //!< individual functions
};

//! Read 1D function from HDF5 dataset
//! \param[in] group HDF5 group containing dataset
//! \param[in] name Name of dataset
//! \return Unique pointer to 1D function
unique_ptr<Function1D> read_function(hid_t group, const char* name);

} // namespace openmc

#endif // OPENMC_ENDF_H
