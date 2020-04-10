#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <mutex>
#include <assert.h>

#include "openmc/openmp_interface.h"

#define NEIGHBOR_SIZE 15

namespace openmc {

//==============================================================================
//! A threadsafe, dynamic container for listing neighboring cells.
//
//! This container is a reduced interface for a linked list with an added OpenMP
//! lock for write operations.  It allows for threadsafe dynamic growth; any
//! number of threads can safely read data without locks or reference counting.
//==============================================================================

class NeighborList
{
public:

  // Attempt to add an element.
  //
  // If the relevant OpenMP lock is currently owned by another thread, this
  // function will return without actually modifying the data.  It has been
  // found that returning the transport calculation and possibly re-adding the
  // element later is slightly faster than waiting on the lock to be released.
  void push_back(int32_t new_elem)
  {
    // Atomically capture the index we want to write to
    int64_t idx;
    #pragma omp atomic capture
    idx = length_++;

    assert(idx < NEIGHBOR_SIZE);

    // Copy element value to the array
    list_[idx] = new_elem;

    #pragma omp atomic write
    status_[idx] = 1;

    // This method has the issue that if some other thread comes in and updates first, length_
    // may be pointing to this functions index, even though it hasn't been written yet.
    /*
    #pragma omp atomic update
    length_++;
    */

    // This method has the issue that if some other thread comes in and updates first, length_
    // may be pointing to one too few. That's ok though, wouldn't result in seg faults, just
    // unnecessary lookups. However, the issue would arise where some other thread was first,
    // but then didn't finish, and this thread now sets length to somewhere past what has been approved...
    /*
    #pragma omp atomic write
    length_ = idx;
    */

    // Could also maintain an array of "ready" variables?
  }

  int64_t get_length() const
  {
    int64_t idx;
    #pragma omp atomic read
    idx = length_;
    assert(idx < NEIGHBOR_SIZE);
    return idx;
  }

  int32_t read_element(int64_t idx) const
  {
    int32_t is_ready;
    #pragma omp atomic read
    is_ready = status_[idx];
    if( is_ready )
      return list_[idx];
    else
      return -1;
  }

private:
  int32_t list_[NEIGHBOR_SIZE];
  int32_t status_[NEIGHBOR_SIZE] = {0};
  int64_t length_ {0};
  //OpenMPMutex mutex_;
  //std::mutex mutex_;
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
