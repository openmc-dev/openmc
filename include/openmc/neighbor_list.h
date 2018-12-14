#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <mutex>

#include "openmc/openmp_interface.h"

namespace openmc {

// Forward declare the neighbor list iterator type.
template <typename T_value, typename T_prefix_iter, typename T_suffix_iter>
class NeighborListIter;

//==============================================================================
//! An efficient, threadsafe, dynamic container for listing neighboring cells.
//
//! This container is a minor improvement upon a linked list.  The linked list
//! allows for threadsafe dynamic growth; any number of threads can safely
//! read from the list without locks or reference counting.  Write access must
//! be protected with a lock, but write events are fairly rare for neighbor
//! lists.  This container also provides the option to flush the linked list
//! into a vector at some convenient time (like between generations where the
//! neighbor lists are not being read) for contiguous data.  The resulting
//! performance benefit is small but consistent.
//==============================================================================

class NeighborList
{
public:
  using value_type = int32_t;
  using prefix_iter = std::vector<value_type>::iterator;
  using suffix_iter = std::forward_list<value_type>::iterator;
  using iterator = NeighborListIter<value_type, prefix_iter, suffix_iter>;

  void push(int new_elem)
  {
    // Try to acquire the lock.
    std::unique_lock<ThreadMutex> lock(mutex_, std::try_to_lock);
    if (lock) {
      // Make sure this element isn't already in the suffix.  Don't check the
      // prefix because we don't expect this function to be called if the "new"
      // element is already there.
      if (std::find(suffix_.cbegin(), suffix_.cend(), new_elem)
          == suffix_.cend()) {
        suffix_.push_front(new_elem);
      }
    }
  }

  void make_consecutive()
  {
    while (!suffix_.empty()) {
      prefix_.push_back(suffix_.front());
      suffix_.pop_front();
    }
  }

  iterator begin();

  iterator end();

private:
  std::vector<value_type> prefix_;
  std::forward_list<value_type> suffix_;
  ThreadMutex mutex_;

  friend class NeighborListIter<value_type, prefix_iter, suffix_iter>;
};

//==============================================================================

template <typename T_value, typename T_prefix_iter, typename T_suffix_iter>
class NeighborListIter
{
public:
  NeighborListIter(NeighborList* nl, T_prefix_iter it)
  {
    base_ = nl;
    if (it != base_->prefix_.end()) {
      in_prefix_ = true;
      prefix_iter_ = it;
    } else {
      in_prefix_ = false;
      suffix_iter_ = base_->suffix_.begin();
    }
  }

  NeighborListIter(NeighborList* nl, T_suffix_iter it)
  {
    in_prefix_ = false;
    suffix_iter_ = it;
  }

  T_value operator*()
  {
    if (in_prefix_) {
      return *prefix_iter_;
    } else {
      return *suffix_iter_;
    }
  }

  bool operator==(const NeighborListIter& other)
  {
    if (in_prefix_ != other.in_prefix_) return false;

    if (in_prefix_) {
      return prefix_iter_ == other.prefix_iter_;
    } else {
      return suffix_iter_ == other.suffix_iter_;
    }
  }

  bool operator!=(const NeighborListIter& other)
  {return !(*this == other);}

  NeighborListIter& operator++()
  {
    if (in_prefix_) {
      // We are in the prefix so increment the prefix iterator.
      ++prefix_iter_;

      // If we've reached the end of the prefix, switch to the suffix iterator.
      if (prefix_iter_ == base_->prefix_.end()) {
        in_prefix_ = false;
        suffix_iter_ = base_->suffix_.begin();
      }

    } else {
      // We are in the suffix so increment the suffix iterator.
      ++suffix_iter_;
    }
    return *this;
  }

private:
  NeighborList* base_;

  union {
    T_prefix_iter prefix_iter_;
    T_suffix_iter suffix_iter_;
  };

  bool in_prefix_;
};

//==============================================================================

inline NeighborList::iterator
NeighborList::begin()
{
  return iterator(this, prefix_.begin());
}

inline NeighborList::iterator
NeighborList::end()
{
  return iterator(this, suffix_.end());
} 

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
