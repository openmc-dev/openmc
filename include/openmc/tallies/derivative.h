#ifndef OPENMC_TALLIES_DERIVATIVE_H
#define OPENMC_TALLIES_DERIVATIVE_H

#include "openmc/particle.h"

#include <unordered_map>
#include <vector>

#include "pugixml.hpp"

//==============================================================================
//! Describes a first-order derivative that can be applied to tallies.
//==============================================================================

namespace openmc {

struct TallyDerivative {
  int id;  //!< User-defined identifier
  int variable;  //!< Independent variable (like temperature)
  int diff_material;  //!< Material this derivative is applied to
  int diff_nuclide;  //!< Nuclide this material is applied to
  double flux_deriv;  //!< Derivative of the current particle's weight

  TallyDerivative() {}
  explicit TallyDerivative(pugi::xml_node node);
};

//==============================================================================
// Non-method functions
//==============================================================================

//! Read tally derivatives from a tallies.xml file
void read_tally_derivatives(pugi::xml_node node);

//! Scale the given score by its logarithmic derivative

void
apply_derivative_to_score(const Particle* p, int i_tally, int i_nuclide,
  double atom_density, int score_bin, double& score);

//! Adjust diff tally flux derivatives for a particle scattering event.
//
//! Note that this subroutine will be called after absorption events in
//! addition to scattering events, but any flux derivatives scored after an
//! absorption will never be tallied.  The paricle will be killed before any
//! further tallies are scored.
//
//! \param p The particle being tracked
void score_collision_derivative(const Particle* p);

//! Adjust diff tally flux derivatives for a particle tracking event.
//
//! \param p The particle being tracked
//! \param distance The distance in [cm] traveled by the particle
void score_track_derivative(const Particle* p, double distance);

//! Set the flux derivatives on differential tallies to zero.
void zero_flux_derivs();

} // namespace openmc

//==============================================================================
// Global variables
//==============================================================================

// Explicit vector template specialization declaration of threadprivate variable
// outside of the openmc namespace for the picky Intel compiler.
extern template class std::vector<openmc::TallyDerivative>;

namespace openmc {

namespace model {
extern std::vector<TallyDerivative> tally_derivs;
#pragma omp threadprivate(tally_derivs)
extern std::unordered_map<int, int> tally_deriv_map;
} // namespace model

// Independent variables
//TODO: convert to enum
constexpr int DIFF_DENSITY {1};
constexpr int DIFF_NUCLIDE_DENSITY {2};
constexpr int DIFF_TEMPERATURE {3};

} // namespace openmc

#endif // OPENMC_TALLIES_DERIVATIVE_H
