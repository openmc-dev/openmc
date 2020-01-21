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

// Different independent variables
enum class DerivativeVariable {
  DENSITY,
  NUCLIDE_DENSITY,
  TEMPERATURE
};

struct TallyDerivative {

  DerivativeVariable variable;  //!< Independent variable (like temperature)
  int id;  //!< User-defined identifier
  int diff_material;  //!< Material this derivative is applied to
  int diff_nuclide;  //!< Nuclide this material is applied to

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
void score_collision_derivative(Particle* p);

//! Adjust diff tally flux derivatives for a particle tracking event.
//
//! \param p The particle being tracked
//! \param distance The distance in [cm] traveled by the particle
void score_track_derivative(Particle* p, double distance);

} // namespace openmc

//==============================================================================
// Global variables
//==============================================================================

namespace openmc {

namespace model {
extern std::vector<TallyDerivative> tally_derivs;
extern std::unordered_map<int, int> tally_deriv_map;
} // namespace model

} // namespace openmc

#endif // OPENMC_TALLIES_DERIVATIVE_H
