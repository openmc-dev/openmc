#ifndef OPENMC_DNP_DRIFT_H
#define OPENMC_DNP_DRIFT_H

#include "openmc/error.h"
#include "openmc/field.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/simulation.h"

namespace openmc {

const int DNP_DRIFT_TRANSPORT_MAX_ITER = 1000;
const double DNP_DRIFT_DISTANCE_MIN = 1.0E-13;

enum class Actions { PLACE_AT_INLET, BLOCK_AT_OUTLET, BLOCK_AT_LOCATION };

//! Adjust position of the fission site.
//!
//! \param[inout] y_n Current position
//! \param[in] y_n_minus_1 Previous position
//! \param[in] t Time
//! \param[in] dt Time step
//! \param[in] decay_time Decay time
void _adjust_position(Position& y_n, const Position& y_n_minus_1, double t,
  double dt, double decay_time);

//! Adjust time.
//!
//! \param[inout] t In: time at position b, out: time at position c (in seconds)
//! \param[in] ta Time at position a (in seconds)
//! \param[in] pa Position a
//! \param[in] pb Position b
//! \param[in] pc Position c (between pa and pb)
void _adjust_time(double& t, double ta, const Position& pa, const Position& pb,
  const Position& pc);

//! Transport of a precursor.
//!
//! \param[inout] p Location of the precursor
//! \param[out] t Final time left before the precursor decays (in seconds)
//! \param[in] decay_time Sampled decay time value (in seconds)
//! \param[in] seed Random number generator seed
bool transport_dnp(SourceSite& site, double decay_time, uint64_t* seed);

} // namespace openmc

#endif // OPENMC_DNP_DRIFT_H
