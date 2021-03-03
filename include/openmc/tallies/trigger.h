#ifndef OPENMC_TALLIES_TRIGGER_H
#define OPENMC_TALLIES_TRIGGER_H

#include <string>

#include "pugixml.hpp"

namespace openmc {

//==============================================================================
// Type definitions
//==============================================================================

enum class TriggerMetric {
  variance, relative_error, standard_deviation, not_active
};

//! Stops the simulation early if a desired tally uncertainty is reached.

struct Trigger {
  TriggerMetric metric;  //!< The type of uncertainty (e.g. std dev) measured
  double threshold;  //!< Uncertainty value below which trigger is satisfied
  int score_index;  //!< Index of the relevant score in the tally's arrays
};

//! Stops the simulation early if a desired k-effective uncertainty is reached.

struct KTrigger
{
  TriggerMetric metric {TriggerMetric::not_active};
  double threshold {0.};
};

//==============================================================================
// Global variable declarations
//==============================================================================

//TODO: consider a different namespace
namespace settings {
  extern KTrigger keff_trigger;
}

//==============================================================================
// Non-memeber functions
//==============================================================================

void check_triggers();

} // namespace openmc
#endif // OPENMC_TALLIES_TRIGGER_H
