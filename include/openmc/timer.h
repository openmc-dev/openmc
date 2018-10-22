#ifndef OPENMC_TIMER_H
#define OPENMC_TIMER_H

#include <chrono>

namespace openmc {

//==============================================================================
//! Class for measuring time elapsed
//==============================================================================

class Timer {
public:
  using clock = std::chrono::high_resolution_clock;

  Timer() {};

  //! Start running the timer
  void start ();

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
  double elapsed_; //!< elasped time in [s]
};

//==============================================================================
// Global variables
//==============================================================================

extern Timer time_bank;
extern Timer time_bank_sample;
extern Timer time_bank_sendrecv;

} // namespace openmc

#endif // OPENMC_TIMER_H
