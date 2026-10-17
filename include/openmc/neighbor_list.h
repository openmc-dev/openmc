#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <mutex>

#include "openmc/openmp_interface.h"

namespace openmc {

//==============================================================================
//! A threadsafe, dynamic container for listing neighboring cells.
//
//! This container is a reduced interface for a growable array with an added
//! OpenMP lock for write operations. It allows for threadsafe dynamic growth;
//! any number of threads can safely read data without locks or reference
//! counting.
//
//! Entries are published through an atomic length: a writer fills the next slot
//! and then releases the new length, and a reader acquires the length before
//! reading any slot. That pairing is what makes the lock-free reads well
//! defined. The first INLINE_CAPACITY entries live in the object itself, so a
//! typical cell never allocates and a search reads the entries from the same
//! cache line as the length. Beyond that the entries are copied into a heap
//! block, which doubles in capacity as needed; a replaced block is retained so
//! that a reader still traversing it stays valid, and since the capacity
//! doubles, the retained blocks hold fewer entries than the live one.
//==============================================================================

class NeighborList {
public:
  using value_type = int32_t;

  //! Number of entries held in the object itself, before any allocation
  static constexpr int INLINE_CAPACITY {8};

  //! A snapshot of the entries that were published when it was taken
  class View {
  public:
    View(const value_type* data, int size) : data_ {data}, size_ {size} {}
    const value_type* begin() const { return data_; }
    const value_type* end() const { return data_ + size_; }
    int size() const { return size_; }

  private:
    const value_type* data_;
    int size_;
  };

  NeighborList() = default;
  NeighborList(const NeighborList&) = delete;
  NeighborList& operator=(const NeighborList&) = delete;

  ~NeighborList()
  {
    void* block = block_.load(std::memory_order_relaxed);
    while (block) {
      void* previous = *static_cast<void**>(block);
      std::free(block);
      block = previous;
    }
  }

  //! Take a snapshot of the entries for lock-free traversal
  //
  //! The length is read first: a reader that observes a length also observes
  //! the storage that held that many entries, since a spill to a larger block
  //! is published before the length that counts the entry which caused it.
  View view() const
  {
    int size = length_.load(std::memory_order_acquire);
    void* block = block_.load(std::memory_order_acquire);
    return View {block ? entries(block) : inline_, size};
  }

  // Attempt to add an element.
  //
  // If the relevant OpenMP lock is currently owned by another thread, this
  // function will return without actually modifying the data.  It has been
  // found that returning to the transport calculation and possibly re-adding
  // the element later is slightly faster than waiting on the lock.
  void push_back(value_type new_elem)
  {
    // Try to acquire the lock.
    std::unique_lock<OpenMPMutex> lock(mutex_, std::try_to_lock);
    if (!lock)
      return;

    // Only one thread writes at a time and taking the lock orders this thread
    // against previous writers, so the entries can be read without ordering.
    int size = length_.load(std::memory_order_relaxed);
    void* block = block_.load(std::memory_order_relaxed);
    value_type* data = block ? entries(block) : inline_;

    // It is possible another thread already added this element to the list
    // while this thread was searching for a cell so make sure the given
    // element isn't a duplicate before adding it.
    if (std::find(data, data + size, new_elem) != data + size)
      return;

    if (size == capacity_) {
      int capacity = 2 * capacity_;
      void* new_block =
        std::malloc(sizeof(void*) + capacity * sizeof(value_type));
      if (!new_block)
        return;

      // Chain the old block so that every block is freed on destruction, and
      // copy the entries a reader may already be traversing.
      *static_cast<void**>(new_block) = block;
      value_type* new_data = entries(new_block);
      std::copy(data, data + size, new_data);

      block_.store(new_block, std::memory_order_release);
      capacity_ = capacity;
      data = new_data;
    }

    // Write the element before publishing it, so a reader that sees the new
    // length is guaranteed to see the element too.
    data[size] = new_elem;
    length_.store(size + 1, std::memory_order_release);
  }

private:
  //! Entries follow the pointer to the previous block at the start of a block
  static value_type* entries(void* block)
  {
    return reinterpret_cast<value_type*>(static_cast<void**>(block) + 1);
  }

  value_type inline_[INLINE_CAPACITY]; //!< entries before any allocation
  std::atomic<void*> block_ {nullptr}; //!< spilled entries, chained to older
  std::atomic<int> length_ {0};        //!< entries published in that storage
  int capacity_ {INLINE_CAPACITY};     //!< slots in that storage; writers only
  OpenMPMutex mutex_;
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
