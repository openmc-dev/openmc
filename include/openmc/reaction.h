//! \file reaction.h
//! Data for an incident neutron reaction

#ifndef OPENMC_REACTION_H
#define OPENMC_REACTION_H

#include <string>

#include "hdf5.h"
#include <gsl/gsl>

#include "openmc/reaction_product.h"
#include "openmc/vector.h"

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
  explicit Reaction(hid_t group, const vector<int>& temperatures);

  //! \brief Calculate reaction rate based on group-wise flux distribution
  //
  //! \param[in] i_temp Temperature index
  //! \param[in] energy Energy group boundaries in [eV]
  //! \param[in] flux Flux in each energy group (not normalized per eV)
  //! \param[in] grid Nuclide energy grid
  //! \return Reaction rate
  double collapse_rate(gsl::index i_temp, gsl::span<const double> energy,
    gsl::span<const double> flux, const vector<double>& grid) const;

  //! Cross section at a single temperature
  struct TemperatureXS {
    int threshold;
    vector<double> value;
  };

  int mt_;                           //!< ENDF MT value
  double q_value_;                   //!< Reaction Q value in [eV]
  bool scatter_in_cm_;               //!< scattering system in center-of-mass?
  bool redundant_;                   //!< redundant reaction?
  vector<TemperatureXS> xs_;         //!< Cross section at each temperature
  vector<ReactionProduct> products_; //!< Reaction products
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
