#ifndef OPENMC_OPENMP_INTERFACE_H
#define OPENMC_OPENMP_INTERFACE_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {

//==============================================================================
//! Accessor functions related to number of threads and thread number
//==============================================================================
inline int num_threads()
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline int thread_num()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

//==============================================================================
//! An object used to prevent concurrent access to a piece of data.
//
//! This type meets the C++ "Lockable" requirements.
//==============================================================================

class OpenMPMutex {
public:
  void init()
  {
#ifdef _OPENMP
    omp_init_lock(&mutex_);
#endif
  }

  OpenMPMutex()
  {
    init();
  }

  ~OpenMPMutex()
  {
#ifdef _OPENMP
    omp_destroy_lock(&mutex_);
#endif
  }

  // omp_lock_t objects cannot be deep copied, they can only be shallow
  // copied. Thus, while shallow copying of an omp_lock_t object is
  // completely valid (provided no race conditions exist), true copying
  // of an OpenMPMutex object is not valid due to the action of the
  // destructor. However, since locks are fungible, we can simply replace
  // copying operations with default construction. This allows storage of
  // OpenMPMutex objects within containers that may need to move/copy them
  // (e.g., std::vector). It is left to the caller to understand that
  // copying of OpenMPMutex does not produce two handles to the same mutex,
  // rather, it produces two different mutexes.

  // Copy constructor
  OpenMPMutex(const OpenMPMutex& other)
  {
    init();
  }

  // Copy assignment operator
  OpenMPMutex& operator=(const OpenMPMutex& other)
  {
    return *this;
  }

  //! Lock the mutex.
  //
  //! This function blocks execution until the lock succeeds.
  void lock()
  {
#ifdef _OPENMP
    omp_set_lock(&mutex_);
#endif
  }

  //! Try to lock the mutex and indicate success.
  //
  //! This function does not block.  It returns immediately and gives false if
  //! the lock is unavailable.
  bool try_lock() noexcept
  {
#ifdef _OPENMP
    return omp_test_lock(&mutex_);
#else
    return true;
#endif
  }

  //! Unlock the mutex.
  void unlock() noexcept
  {
#ifdef _OPENMP
    omp_unset_lock(&mutex_);
#endif
  }

private:
#ifdef _OPENMP
  omp_lock_t mutex_;
#endif
};

} // namespace openmc
#endif // OPENMC_OPENMP_INTERFACE_H
