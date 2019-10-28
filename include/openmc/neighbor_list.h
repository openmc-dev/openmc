#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <mutex>

#include "openmc/openmp_interface.h"

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
  using value_type = int32_t;
  using const_iterator = std::forward_list<value_type>::const_iterator;

  // Attempt to add an element.
  //
  // If the relevant OpenMP lock is currently owned by another thread, this
  // function will return without actually modifying the data.  It has been
  // found that returning the transport calculation and possibly re-adding the
  // element later is slightly faster than waiting on the lock to be released.
  void push_back(int new_elem)
  {
    // Try to acquire the lock.
    std::unique_lock<OpenMPMutex> lock(mutex_, std::try_to_lock);
    if (lock) {
      // It is possible another thread already added this element to the list
      // while this thread was searching for a cell so make sure the given
      // element isn't a duplicate before adding it.
      if (std::find(list_.cbegin(), list_.cend(), new_elem) == list_.cend()) {
        // Find the end of the list and add the the new element there.
        if (!list_.empty()) {
          auto it1 = list_.cbegin();
          auto it2 = ++list_.cbegin();
          while (it2 != list_.cend()) it1 = it2++;
          list_.insert_after(it1, new_elem);
        } else {
          list_.push_front(new_elem);
        }
      }
    }
  }

  const_iterator cbegin() const
  {return list_.cbegin();}

  const_iterator cend() const
  {return list_.cend();}


private:
  std::forward_list<value_type> list_;
  OpenMPMutex mutex_;
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
