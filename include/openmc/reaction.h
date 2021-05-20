//! \file reaction.h
//! Data for an incident neutron reaction

#ifndef OPENMC_REACTION_H
#define OPENMC_REACTION_H

#include <string>
#include <vector>

#include <gsl/gsl>
#include "hdf5.h"

#include "openmc/particle.h"
#include "openmc/reaction_product.h"
#include "openmc/serialize.h"

namespace openmc {

//==============================================================================
//! Data for a single reaction including cross sections (possibly at multiple
//! temperatures) and reaction products (with secondary angle-energy
//! distributions)
//==============================================================================

class Reaction {
public:
  //! Construct reaction from HDF5 data
  //! \param[in] group HDF5 group containing reaction data
  //! \param[in] temperatures Desired temperatures for cross sections
  explicit Reaction(hid_t group, const std::vector<int>& temperatures);

  double xs(gsl::index i_temp, gsl::index i_grid, double interp_factor) const;
  double xs(const NuclideMicroXS& micro) const;

  //! \brief Calculate reaction rate based on group-wise flux distribution
  //
  //! \param[in] i_temp Temperature index
  //! \param[in] energy Energy group boundaries in [eV]
  //! \param[in] flux Flux in each energy group (not normalized per eV)
  //! \param[in] grid Nuclide energy grid
  //! \return Reaction rate
  double collapse_rate(gsl::index i_temp, gsl::span<const double> energy,
    gsl::span<const double> flux, const std::vector<double>& grid) const;

  void serialize(DataBuffer& buffer) const;

  //! Cross section at a single temperature
  struct TemperatureXS {
    int threshold;
    std::vector<double> value;
  };

  int mt_;             //!< ENDF MT value
  double q_value_;     //!< Reaction Q value in [eV]
  bool scatter_in_cm_; //!< scattering system in center-of-mass?
  bool redundant_;     //!< redundant reaction?
  std::vector<TemperatureXS> xs_; //!< Cross section at each temperature
  std::vector<ReactionProduct> products_; //!< Reaction products
};

class ReactionFlat {
public:
  // Constructors
  #pragma omp declare target
  explicit ReactionFlat(const uint8_t* data);
  #pragma omp end declare target

  double xs(gsl::index i_temp, gsl::index i_grid, double interp_factor) const;
  double xs(const NuclideMicroXS& micro) const;

  double collapse_rate(gsl::index i_temp, gsl::span<const double> energy,
    gsl::span<const double> flux, const std::vector<double>& grid) const;

  // Accessors
  int mt() const;
  double q_value() const;
  bool scatter_in_cm() const;
  bool redundant() const;
  int xs_threshold(gsl::index i_temp) const;
  gsl::span<const double> xs_value(gsl::index i_temp) const;
  ReactionProductFlat products(gsl::index i) const;
  size_t n_xs() const { return n_xs_; }
  size_t n_products() const { return n_products_; }
private:
  const uint8_t* data_;
  size_t n_xs_;
  size_t n_products_;
};

class ReactionFlatContainer {
public:
  // Constructors
  explicit ReactionFlatContainer(const Reaction& rx);

  void copy_to_device();
  void release_from_device();

  ReactionFlat obj() const;
  int mt() const { return this->obj().mt(); }
private:
  // Data members
  DataBuffer buffer_;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Return reaction name given an ENDF MT value
//
//! \param[in] mt  ENDF MT value
//! \return Name of the corresponding reaction
std::string reaction_name(int mt);

//! Return reaction type (MT value) given a reaction name
//
//! \param[in] name  Reaction name
//! \return Corresponding reaction type (MT value)
int reaction_type(std::string name);

} // namespace openmc

#endif // OPENMC_REACTION_H
