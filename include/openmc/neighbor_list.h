#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <mutex>
#include <assert.h>

#include "openmc/openmp_interface.h"

#define NEIGHBOR_SIZE 11 // limited by triso regression test

namespace openmc{

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

  NeighborList()
  {
    for(int i = 0; i < NEIGHBOR_SIZE; i++)
      list_[i] = -1;
  }

  // Attempt to add an element.
  //
  // If the relevant OpenMP lock is currently owned by another thread, this
  // function will return without actually modifying the data.  It has been
  // found that returning the transport calculation and possibly re-adding the
  // element later is slightly faster than waiting on the lock to be released.
  void push_back(int32_t new_elem)
  {
    printf("pushing back %d\n", new_elem);
    for( int i = 0; i < NEIGHBOR_SIZE; i++)
    {
      // This line checks to see if the new_elem is already in the list
      int retrieved_id;

      // OpenMP 5.1
      /*
         #pragma omp atomic compare capture
         {
         retrieved_id = list_[i];
         if (list_[i] == -1)
         list_[i] = new_elem;
         }
      */

      // Compiler non-portable builtin atomic CAS. Unclear if this is sequentially consistent or not (needs to be!)
      __sync_synchronize();
      retrieved_id = __sync_val_compare_and_swap(&list_[i], -1, new_elem);
      __sync_synchronize();

      // Case 1: The element was not initialized yet, so the previous line had the effect of setting it to new_elem and returning -1.
      // Case 2: The element was already initialized to the current new_elem, so the atomicCAS call will return new_elem
      if( retrieved_id == -1 || retrieved_id == new_elem)
        return;

      // Case 3: The element was already initialized to a different cell_id, so it will return some other value != -1 and != new_elem
      // so, we continue reading through the list.
    }
  }
  
  int32_t list_[NEIGHBOR_SIZE];
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
