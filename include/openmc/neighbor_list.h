#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <iostream>
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

// First, the implementation for GPU/CPU execution which uses a little more
// memory than the host version, and secondly, the original host-only version.
// Linked lists are very poorly performing on GPU, so the GPU-compatible
// uses a flagged buffer approach. If the buffer runs out of space, the kernel
// simply stops appending to the list until its done. After running the
// geometry kernel, the host checks if any buffers had attempts to write
// past their end. The buffers are then expanded on the host as necessary.
//
// In addition, it's possible for duplicate neighbors to be written since
// a general locking method is really not possible in CUDA. As a result,
// another kernel checks for duplicates after completion of this kernel.
//
// The neighbor lists thus fill out a bit more slowly than in the pure
// CPU version, but asymptotically, the performance is the same.

#ifdef __CUDACC__
class NeighborList {
public:
  using value_type = int32_t;
  using const_iterator = int32_t const*;
  using iterator = int32_t*;
  static constexpr unsigned buffer_size = 11;

  __host__ __device__ NeighborList() : lock(0), current_size(0) {}

  __host__ __device__ void push_back(int new_elem)
  {
    // Ensure that the new_elem is not already in the list. This is pretty
    // divergent code, but it really only needs to run a few times at the
    // very start of a calculation.
#ifdef __CUDA_ARCH__
    // Locking below acts as a mutex on blocks. The loop
    // here acts as a mutex on threads. Both HAVE to be
    // treated separately in CUDA.
    for (unsigned thread_id = 0; thread_id < blockDim.x; ++thread_id) {

      if (threadIdx.x == thread_id) {

        // Blockwise lock
        while (atomicCAS(&lock, 0, 1) != 0)
          ;

        // Check if this thing we're pushing back was put in by someone else
        for (unsigned i = 0; i < current_size; ++i) {
          if (new_elem == buffer_[i]) {
            atomicExch(&lock, 0);
            return;
          }
        }

        // Fail if writing out of buffer
        if (current_size == buffer_size)
          asm("trap;");

        // Add in the new element
        buffer_[current_size++] = new_elem;

        // Unlock block lock
        atomicExch(&lock, 0);
      }
    }
#else
// Ahhh, so nice to be able to use critical sections on the CPU
#pragma omp critical
    {
      if (buffer_size == current_size) {
        std::cerr << "attempt to write past end of neighbor list buffer."
                  << std::endl
                  << "Must recompile with increased buffer size." << std::endl;
      }
      buffer_[current_size++] = new_elem;
    }
#endif
  }

  __host__ __device__ const_iterator cbegin() const { return buffer_; }
  __host__ __device__ const_iterator cend() const
  {
    return buffer_ + current_size;
  }
  __host__ __device__ iterator begin() { return buffer_; }
  __host__ __device__ iterator end() { return buffer_ + current_size; }
  __host__ __device__ const_iterator begin() const { return buffer_; }
  __host__ __device__ const_iterator end() const
  {
    return buffer_ + current_size;
  }

private:
  int lock;
  unsigned int current_size;
  value_type buffer_[buffer_size]; // _temporarily_ hardcoded
};

#else
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
#endif

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
