#include "openmc/random_ray/tally_convert.h"
#include "openmc/geometry.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace random_ray {

std::vector<std::vector<TallyTask>> tally_task;

}

//==============================================================================
// Non-method functions
//==============================================================================


} // namespace openmc
