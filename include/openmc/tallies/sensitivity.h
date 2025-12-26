#ifndef OPENMC_TALLIES_SENSITIVITY_H
#define OPENMC_TALLIES_SENSITIVITY_H

#include "openmc/particle.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/vector.h"

#include <unordered_map>

#include "pugixml.hpp"

//==============================================================================
//! Describes a first-order sensitivity that can be applied to importance 
//  weighted tallies.
//==============================================================================

namespace openmc {

//==============================================================================
// Class Definitions
//==============================================================================

class SensitivityTally : public Tally
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions
  SensitivityTally(int32_t id);
  SensitivityTally(pugi::xml_node node);
  virtual ~SensitivityTally();
  //static SensitivityTally* create(int32_t id = -1);

  //----------------------------------------------------------------------------
  // Accessors

  void set_filters(span<Filter*> filters) override;

  //----------------------------------------------------------------------------
  // Other methods.
  void add_filter(Filter* filter) override { set_filters({&filter, 1}); }
  
  void init_results() override;

  void accumulate() override;

  void reset() override;

  //----------------------------------------------------------------------------
  // Major public data members.

  double denominator_ {0.0}; //!< Denominator of the sensitivity, <\phi F \psi>

  xt::xtensor<double, 3> previous_results_;

  bool gpt_ {false}; //!< if true, do not use denominator.

  //int sens_ {C_NONE}; //!< Index of a Sensitivity object for sensitivity tallies.
};

// Different independent variables
enum class SensitivityVariable {
  CROSS_SECTION,
  MULTIPOLE,
  CURVE_FIT
};

struct TallySensitivity {

  SensitivityVariable variable;  //!< Independent variable (like xs)
  int id;  //!< User-defined identifier
  int sens_nuclide;     //!< Nuclide the sensitivity is calculated for
  int sens_element;     //!< Element the sensitivity is calculated for
  int sens_reaction;    //!< Need something to specify reaction, use ReactionType?
  int sens_material;    //!< Material the sensitivity is calculated in
  int sens_cell;        //!< Cell the sensitivity is calculated in
  std::vector<double> energy_bins_; //!< Energy bins on which to discretize the cross section
  int n_bins_; //!< something to indicate the size of the vector

  TallySensitivity() {}
  explicit TallySensitivity(pugi::xml_node node);

  void set_bins(span<const double> bins);
};

//==============================================================================
// Non-method functions
//==============================================================================

//! Read tally sensitivities from a tallies.xml file
void read_tally_sensitivities(pugi::xml_node node);

//! Scale the given score by its logarithmic derivative

//void
//apply_sensitivity_to_score(const Particle* p, int i_tally, int i_nuclide,
//  double atom_density, int score_bin, double& score);

//! Adjust diff tally flux derivatives for a particle scattering event. 
//
//! Note that this subroutine will be called after absorption events in
//! addition to scattering events, but any flux derivatives scored after an
//! absorption will never be tallied.  The paricle will be killed before any
//! further tallies are scored.
//
//! \param p The particle being tracked
void score_collision_sensitivity(Particle& p);

//! Adjust diff tally flux derivatives for a particle tracking event.
//
//! \param p The particle being tracked
//! \param distance The distance in [cm] traveled by the particle
void score_track_sensitivity(Particle& p, double distance);

void score_source_sensitivity(Particle& p);

} // namespace openmc

//==============================================================================
// Global variables
//==============================================================================

namespace openmc {

namespace model {
extern std::vector<TallySensitivity> tally_sens;
extern std::unordered_map<int, int> tally_sens_map;
} // namespace model

} // namespace openmc

#endif // OPENMC_TALLIES_SENSITIVITY_H