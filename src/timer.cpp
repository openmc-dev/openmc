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

Timer time_bank;
Timer time_bank_sample;
Timer time_bank_sendrecv;

//==============================================================================
// Fortran compatibility
//==============================================================================

extern "C" double time_bank_elapsed() { return time_bank.elapsed(); }
extern "C" double time_bank_sample_elapsed() { return time_bank_sample.elapsed(); }
extern "C" double time_bank_sendrecv_elapsed() { return time_bank_sendrecv.elapsed(); }

extern "C" void reset_timers()
{
  time_bank.reset();
  time_bank_sample.reset();
  time_bank_sendrecv.reset();
}

} // namespace openmc
