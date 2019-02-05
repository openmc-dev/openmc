#include "openmc/timer.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

Timer time_active;
Timer time_bank;
Timer time_bank_sample;
Timer time_bank_sendrecv;
Timer time_finalize;
Timer time_inactive;
Timer time_initialize;
Timer time_read_xs;
Timer time_tallies;
Timer time_total;
Timer time_transport;

} // namespace simulation

//==============================================================================
// Timer implementation
//==============================================================================

void Timer::start ()
{
  running_ = true;
  start_ = clock::now();
}

void Timer::stop()
{
  elapsed_ = elapsed();
  running_ = false;
}

void Timer::reset()
{
  running_ = false;
  elapsed_ = 0.0;
}

double Timer::elapsed()
{
  if (running_) {
    std::chrono::duration<double> diff = clock::now() - start_;
    return elapsed_ + diff.count();
  } else {
    return elapsed_;
  }
}

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" double time_active_elapsed() { return simulation::time_active.elapsed(); }
extern "C" double time_bank_elapsed() { return simulation::time_bank.elapsed(); }
extern "C" double time_bank_sample_elapsed() { return simulation::time_bank_sample.elapsed(); }
extern "C" double time_bank_sendrecv_elapsed() { return simulation::time_bank_sendrecv.elapsed(); }
extern "C" double time_finalize_elapsed() { return simulation::time_finalize.elapsed(); }
extern "C" double time_inactive_elapsed() { return simulation::time_inactive.elapsed(); }
extern "C" double time_initialize_elapsed() { return simulation::time_initialize.elapsed(); }
extern "C" double time_read_xs_elapsed() { return simulation::time_read_xs.elapsed(); }
extern "C" double time_tallies_elapsed() { return simulation::time_tallies.elapsed(); }
extern "C" double time_total_elapsed() { return simulation::time_total.elapsed(); }
extern "C" double time_transport_elapsed() { return simulation::time_transport.elapsed(); }

//==============================================================================
// Non-member functions
//==============================================================================

void reset_timers()
{
  simulation::time_active.reset();
  simulation::time_bank.reset();
  simulation::time_bank_sample.reset();
  simulation::time_bank_sendrecv.reset();
  simulation::time_finalize.reset();
  simulation::time_inactive.reset();
  simulation::time_initialize.reset();
  simulation::time_read_xs.reset();
  simulation::time_tallies.reset();
  simulation::time_total.reset();
  simulation::time_transport.reset();
}

} // namespace openmc
