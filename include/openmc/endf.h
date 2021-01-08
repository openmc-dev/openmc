//! \file endf.h
//! Classes and functions related to the ENDF-6 format

#ifndef OPENMC_ENDF_H
#define OPENMC_ENDF_H

#include <memory>
#include <vector>

#include "hdf5.h"
#include <gsl/gsl>

#include "openmc/constants.h"
#include "openmc/serialize.h"

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
  virtual void serialize(DataBuffer& buffer) const = 0;
};

//==============================================================================
//! One-dimensional function expressed as a polynomial
//==============================================================================

class Polynomial : public Function1D {
public:
  //! Construct polynomial from HDF5 data
  //! \param[in] dset Dataset containing coefficients
  explicit Polynomial(hid_t dset);

  //! Evaluate the polynomials
  //! \param[in] x independent variable
  //! \return Polynomial evaluated at x
  double operator()(double x) const override;

  void serialize(DataBuffer& buffer) const;
private:
  std::vector<double> coef_; //!< Polynomial coefficients
};

class PolynomialFlat {
public:
  explicit PolynomialFlat(const uint8_t* data) : data_(data) { }

  double operator()(double x) const;
private:
  gsl::span<const double> coef() const;

  const uint8_t* data_;
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

  //! Evaluate the tabulated function
  //! \param[in] x independent variable
  //! \return Function evaluated at x
  double operator()(double x) const override;

  void serialize(DataBuffer& buffer) const;

  // Accessors
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& y() const { return y_; }
private:
  std::size_t n_regions_ {0}; //!< number of interpolation regions
  std::vector<int> nbt_; //!< values separating interpolation regions
  std::vector<Interpolation> int_; //!< interpolation schemes
  std::size_t n_pairs_; //!< number of (x,y) pairs
  std::vector<double> x_; //!< values of abscissa
  std::vector<double> y_; //!< values of ordinate
};

class Tabulated1DFlat {
public:
  explicit Tabulated1DFlat(const uint8_t* data);

  double operator()(double x) const;

private:
  gsl::span<const int> nbt() const;
  Interpolation interp(gsl::index i) const;
  gsl::span<const double> x() const;
  gsl::span<const double> y() const;

  const uint8_t* data_;
  size_t n_regions_;
  size_t n_pairs_;
};

//==============================================================================
//! Coherent elastic scattering data from a crystalline material
//==============================================================================

class CoherentElasticXS : public Function1D {
public:
  explicit CoherentElasticXS(hid_t dset);

  double operator()(double E) const override;

  const std::vector<double>& bragg_edges() const { return bragg_edges_; }
  const std::vector<double>& factors() const { return factors_; }

  void serialize(DataBuffer& buffer) const;
private:
  std::vector<double> bragg_edges_; //!< Bragg edges in [eV]
  std::vector<double> factors_;     //!< Partial sums of structure factors [eV-b]
};

class CoherentElasticXSFlat {
public:
  explicit CoherentElasticXSFlat(const uint8_t* data) : data_{data} { }

  double operator()(double E) const;

  gsl::span<const double> bragg_edges() const;
  gsl::span<const double> factors() const;
private:
  const uint8_t* data_;
};

//==============================================================================
//! Incoherent elastic scattering cross section
//==============================================================================

class IncoherentElasticXS : public Function1D {
public:
  explicit IncoherentElasticXS(hid_t dset);

  double operator()(double E) const override;

  void serialize(DataBuffer& buffer) const;
private:
  double bound_xs_; //!< Characteristic bound xs in [b]
  double debye_waller_; //!< Debye-Waller integral divided by atomic mass in [eV^-1]
};

class IncoherentElasticXSFlat {
public:
  explicit IncoherentElasticXSFlat(const uint8_t* data) : data_{data} { }

  double operator()(double E) const;
private:
  double bound_xs() const { return *reinterpret_cast<const double*>(data_ + 4); }
  double debye_waller() const { return *reinterpret_cast<const double*>(data_ + 12); }

  const uint8_t* data_;
};

//! Read 1D function from HDF5 dataset
//! \param[in] group HDF5 group containing dataset
//! \param[in] name Name of dataset
//! \return Unique pointer to 1D function
std::unique_ptr<Function1D> read_function(hid_t group, const char* name);

} // namespace openmc

#endif // OPENMC_ENDF_H
