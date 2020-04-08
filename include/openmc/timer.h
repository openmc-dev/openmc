#ifndef OPENMC_TIMER_H
#define OPENMC_TIMER_H

#include <chrono>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Timer;

namespace simulation {

extern Timer time_active;
extern Timer time_bank;
extern Timer time_bank_sample;
extern Timer time_bank_sendrecv;
extern Timer time_finalize;
extern Timer time_inactive;
extern Timer time_initialize;
extern Timer time_read_xs;
extern Timer time_sample_source;
extern Timer time_tallies;
extern Timer time_total;
extern Timer time_transport;
extern Timer time_event_init;
extern Timer time_event_calculate_xs;
extern Timer time_event_advance_particle;
extern Timer time_event_surface_crossing;
extern Timer time_event_collision;
extern Timer time_event_death;

} // namespace simulation

//==============================================================================
//! Class for measuring time elapsed
//==============================================================================

class Timer {
public:
  using clock = std::chrono::high_resolution_clock;

  Timer() {};

  //! Start running the timer
  void start();

  //! Get total elapsed time in seconds
  //! \return Elapsed time in [s]
  double elapsed();

  //! Stop running the timer
  void stop();

  //! Stop the timer and reset its elapsed time
  void reset();

private:
  bool running_ {false}; //!< is timer running?
  std::chrono::time_point<clock> start_; //!< starting point for clock
  double elapsed_ {0.0}; //!< elapsed time in [s]
};

//==============================================================================
// Non-member functions
//==============================================================================

void reset_timers();

} // namespace openmc

#endif // OPENMC_TIMER_H
