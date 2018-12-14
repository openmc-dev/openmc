#ifndef OPENMC_OPENMP_INTERFACE_H
#define OPENMC_OPENMP_INTERFACE_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {

class ThreadMutex
{

#ifdef _OPENMP
//==============================================================================
// Implementation for when OpenMP is actually defined.

public:
  ThreadMutex()
  {omp_init_lock(&mutex_);}

  ~ThreadMutex()
  {omp_destroy_lock(&mutex_);}

  void lock()
  {omp_set_lock(&mutex_);}

  bool try_lock()
  {return omp_test_lock(&mutex_);}

  void unlock()
  {omp_unset_lock(&mutex_);}

private:
  omp_lock_t mutex_;

#else
//==============================================================================
// Empty implementation for when OpenMP is not defined.

public:
  ThreadMutex() {}
  ~ThreadMutex() = default;
  void lock() {}
  bool try_lock() {return true;}
  void unlock() {}
#endif
};

} // namespace openmc
#endif // OPENMC_OPENMP_INTERFACE_H
