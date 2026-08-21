#ifndef OPENMC_DNP_DRIFT_H
#define OPENMC_DNP_DRIFT_H

#include "openmc/particle_data.h"
#include "openmc/position.h"

namespace openmc {

const int DNP_DRIFT_TRANSPORT_MAX_ITER =
  1000; //! Maximum number of iterations allowed in the DNP transport loop.
const double DNP_DRIFT_DISTANCE_MIN =
  1.0E-13; //! Minimum distance between two integration steps before stopping
           //! the DNP transport loop.

enum class Actions { REENTER, ESCAPE, DECAY_IN_PLACE };

//! Adjust position linearly.
//!
//! This function is used to adjust the position of a DNP that reached its decay
//! time between two integration steps (i.e., t is greater than decay_time). The
//! final position of the DNP is adjusted linearly so that the time
//! corresponding to this new position is equal to the sampled decay time.
//!
//! \param[inout] y_n Current position
//! \param[in] y_n_minus_1 Previous position
//! \param[in] t Time at current position
//! \param[in] dt Time step
//! \param[in] decay_time Decay time
void _adjust_position(Position& y_n, const Position& y_n_minus_1, double t,
  double dt, double decay_time);

//! Adjust time linearly.
//!
//! This function is used to adjust the time associated with the position of a
//! DNP that needs to be stopped between two integration steps. The position C,
//! where the DNP is stopped, must be located between position A and B. Time is
//! adjusted linearly.
//!
//! \param[inout] t In: time at position B, out: time at position C
//! \param[in] ta Time at position A
//! \param[in] pa Position A
//! \param[in] pb Position B
//! \param[in] pc Position C (between A and B)
void _adjust_time(double& t, double ta, const Position& pa, const Position& pb,
  const Position& pc);

//! Explicitly transport a delayed neutron precursor using streamline
//! integration.
//!
//! \param[inout] site Fission site corresponding to the DNP
//! \param[in] decay_time Sampled decay time (in seconds)
//! \param[in] seed Random number generator seed
//! \return True if the site is considered inside the explicitly modeled part of
//!         the system. False otherwise.
bool transport_dnp(SourceSite& site, double decay_time, uint64_t* seed);

} // namespace openmc

#endif // OPENMC_DNP_DRIFT_H
