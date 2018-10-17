#include "openmc/timer.h"

namespace openmc {

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
    return elapsed_ += diff.count();
  } else {
    return elapsed_;
  }
}

//==============================================================================
// Global variables
//==============================================================================

Timer time_active;
Timer time_bank;
Timer time_bank_sample;
Timer time_bank_sendrecv;
Timer time_finalize;
Timer time_inactive;
Timer time_initialize;
Timer time_total;

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" double time_active_elapsed() { return time_active.elapsed(); }
extern "C" double time_bank_elapsed() { return time_bank.elapsed(); }
extern "C" double time_bank_sample_elapsed() { return time_bank_sample.elapsed(); }
extern "C" double time_bank_sendrecv_elapsed() { return time_bank_sendrecv.elapsed(); }
extern "C" double time_finalize_elapsed() { return time_finalize.elapsed(); }
extern "C" double time_inactive_elapsed() { return time_inactive.elapsed(); }
extern "C" double time_initialize_elapsed() { return time_initialize.elapsed(); }
extern "C" double time_total_elapsed() { return time_total.elapsed(); }

extern "C" void reset_timers()
{
  time_active.reset();
  time_bank.reset();
  time_bank_sample.reset();
  time_bank_sendrecv.reset();
  time_finalize.reset();
  time_initialize.reset();
  time_total.reset();
}

} // namespace openmc
