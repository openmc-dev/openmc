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
  HD virtual xsfloat operator()(xsfloat x) const = 0;
  virtual ~Function1D() = default;
  Function1D() = default;
  Function1D(Function1D const&) = delete;
  Function1D(Function1D&&) = default;
  Function1D& operator=(Function1D&&) = default;
};

//==============================================================================
//! One-dimensional function expressed as a polynomial
//==============================================================================

class Polynomial : public Function1D {
public:
  //! Construct polynomial from HDF5 data
  //! \param[in] dset Dataset containing coefficients
  explicit Polynomial(hid_t dset);
  Polynomial(Polynomial&&) = default;

  //! Evaluate the polynomials
  //! \param[in] x independent variable
  //! \return Polynomial evaluated at x
  HD xsfloat operator()(xsfloat x) const override;
private:
  vector<xsfloat> coef_; //!< Polynomial coefficients
};

//==============================================================================
//! One-dimensional interpolable function
//==============================================================================

class Tabulated1D : public Function1D {
public:
  Tabulated1D() = default;
  Tabulated1D& operator=(Tabulated1D&&) = default;
  Tabulated1D(Tabulated1D&&) = default;

  //! Construct function from HDF5 data
  //! \param[in] dset Dataset containing tabulated data
  explicit Tabulated1D(hid_t dset);

  //! Evaluate the tabulated function
  //! \param[in] x independent variable
  //! \return Function evaluated at x
  HD xsfloat operator()(xsfloat x) const override;

  // Accessors
  const vector<xsfloat>& x() const { return x_; }
  const vector<xsfloat>& y() const { return y_; }

private:
  std::size_t n_regions_ {0}; //!< number of interpolation regions
  vector<int> nbt_;           //!< values separating interpolation regions
  vector<Interpolation> int_; //!< interpolation schemes
  std::size_t n_pairs_; //!< number of (x,y) pairs
  vector<xsfloat> x_;    //!< values of abscissa
  vector<xsfloat> y_;    //!< values of ordinate
};

//==============================================================================
//! Coherent elastic scattering data from a crystalline material
//==============================================================================

class CoherentElasticXS : public Function1D {
public:
  explicit CoherentElasticXS(hid_t dset);
  CoherentElasticXS(CoherentElasticXS&&) = default;

  HD xsfloat operator()(xsfloat E) const override;

  HD const vector<xsfloat>& bragg_edges() const { return bragg_edges_; }
  HD const vector<xsfloat>& factors() const { return factors_; }

private:
  vector<xsfloat> bragg_edges_; //!< Bragg edges in [eV]
  vector<xsfloat> factors_;     //!< Partial sums of structure factors [eV-b]
};

//==============================================================================
//! Incoherent elastic scattering cross section
//==============================================================================

class IncoherentElasticXS : public Function1D {
public:
  explicit IncoherentElasticXS(hid_t dset);
  IncoherentElasticXS(IncoherentElasticXS&&) = default;

  HD xsfloat operator()(xsfloat E) const override;
private:
  double bound_xs_; //!< Characteristic bound xs in [b]
  double debye_waller_; //!< Debye-Waller integral divided by atomic mass in [eV^-1]
};

//! Read 1D function from HDF5 dataset
//! \param[in] group HDF5 group containing dataset
//! \param[in] name Name of dataset
//! \return Unique pointer to 1D function
unique_ptr<Function1D> read_function(hid_t group, const char* name);

} // namespace openmc

#endif // OPENMC_ENDF_H
