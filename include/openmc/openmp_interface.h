#ifndef OPENMC_OPENMP_INTERFACE_H
#define OPENMC_OPENMP_INTERFACE_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {

//==============================================================================
//! An object used to prevent concurrent access to a piece of data.
//
//! This type meets the C++ "Lockable" requirements.
//==============================================================================

class OpenMPMutex
{
public:
  OpenMPMutex()
  {
    #ifdef _OPEMP
      omp_init_lock(&mutex_);
    #endif
  }

  ~OpenMPMutex()
  {
    #ifdef _OPEMP
      omp_destroy_lock(&mutex_);
    #endif
  }

  // Mutexes cannot be copied.  We need to explicitly delete the copy
  // constructor and copy assignment operator to ensure the compiler doesn't
  // "help" us by implicitly trying to copy the underlying mutexes.
  OpenMPMutex(const OpenMPMutex&) = delete;
  OpenMPMutex& operator= (const OpenMPMutex&) = delete;

  //! Lock the mutex.
  //
  //! This function blocks execution until the lock succeeds.
  void lock()
  {
    #ifdef _OPEMP
      omp_set_lock(&mutex_);
    #endif
  }

  //! Try to lock the mutex and indicate success.
  //
  //! This function does not block.  It returns immediately and gives false if
  //! the lock is unavailable.
  bool try_lock() noexcept
  {
    #ifdef _OPEMP
      return omp_test_lock(&mutex_);
    #else
      return true;
    #endif
  }

  //! Unlock the mutex.
  void unlock() noexcept
  {
    #ifdef _OPEMP
      omp_unset_lock(&mutex_);
    #endif
  }

private:
  #ifdef _OPEMP
    omp_lock_t mutex_;
  #endif
};

} // namespace openmc
#endif // OPENMC_OPENMP_INTERFACE_H
