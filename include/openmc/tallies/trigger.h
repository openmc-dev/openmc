#ifndef OPENMC_TALLIES_TRIGGER_H
#define OPENMC_TALLIES_TRIGGER_H

#include <string>

#include "pugixml.hpp"

namespace openmc {

//! Stops the simulation early if a desired tally uncertainty is reached.

extern "C" struct Trigger
{
  int type;  //!< variance, std_dev, or rel_err
  double threshold;  //!< uncertainty value below which trigger is satisfied
  int score_index;  //!< index of the relevant score in the tally's arrays
  double variance {0.};
  double std_dev {0.};
  double rel_err {0.};
};

//! Stops the simulation early if a desired k-effective uncertainty is reached.

extern "C" struct KTrigger
{
  int type{0};
  double threshold {0.};
};

//TODO: consider a different namespace
namespace settings {
  extern "C" KTrigger keff_trigger;
}

} // namespace openmc
#endif // OPENMC_TALLIES_TRIGGER_H
