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
    // Lock the object
    mutex_.lock();

    // Check to see if the element already exists in the list
    for( int64_t i = 0; i < length_; i++ )
    {
      if( list_[i] == new_elem )
      {
        mutex_.unlock();
        return;
      }
    }

    // Determine the index we want to write to
    int64_t idx = length_;
    assert(idx < NEIGHBOR_SIZE);
    
    // Copy element value to the array
    list_[idx] = new_elem;

    #pragma omp atomic
    length_++;

    // unlock the object
    mutex_.unlock();
  }

  int64_t get_length() const
  {
    int64_t idx;
    #pragma omp atomic read
    idx = length_;
    assert(idx < NEIGHBOR_SIZE);
    return idx;
  }
  
  int32_t list_[NEIGHBOR_SIZE];

private:
  int64_t length_ {0};
  OpenMPMutex mutex_;
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
