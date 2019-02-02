#ifndef OPENMC_TALLIES_TRIGGER_H
#define OPENMC_TALLIES_TRIGGER_H

#include <string>

#include "pugixml.hpp"

namespace openmc {

//==============================================================================
// Structs
//==============================================================================

//! Stops the simulation early if a desired tally uncertainty is reached.

struct Trigger
{
  Trigger(int type_, double threshold_, int score_index_)
    : type(type_), threshold(threshold_), score_index(score_index_) {}

  int type;  //!< variance, std_dev, or rel_err
  double threshold;  //!< uncertainty value below which trigger is satisfied
  int score_index;  //!< index of the relevant score in the tally's arrays
  double variance {0.};
  double std_dev {0.};
  double rel_err {0.};
};

//! Stops the simulation early if a desired k-effective uncertainty is reached.

struct KTrigger
{
  int type {0};
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
