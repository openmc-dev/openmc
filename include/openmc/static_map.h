#ifndef OPENMC_STATIC_MAP_H
#define OPENMC_STATIC_MAP_H

#include <algorithm> // for stable_sort, find_if, lower_bound
#include <utility> // for pair

#include "openmc/vector.h"

namespace openmc {

//==========================================================================
// Integer hash function
//==========================================================================

//! Hash function for integers that simply returns the integer.
//
//! If the values in the hash table are consecutive starting from zero, this
//! will produce no collisions. Otherwise, all bets are off.
//
//! TODO: explore the performance of different hash functions
template<typename Key>
struct integer_hash
{
  using value_type = Key;
  using result_type = std::size_t;

  result_type operator()(value_type value) const
  {
    return static_cast<result_type>(value);
  }
};

//==========================================================================
// Class declarations
//==========================================================================

//! Fixed-size map with on-device access.
//
//! This is a read-only hash table that supports on-device lookups. The table
//! is constructed on the host and must be "finalized" before it can be used to
//! search for items. Adapted from the \c Static_Device_Map class in Oak Ridge
//! National Laboratory's SCALE code system with help from 
//! Seth Johnson, Tom Evans, and Steven Hamilton.
template<typename Key, typename T, typename Hash = integer_hash<Key>>
class static_map {
public:
  // Types, aliases
  using key_type = Key;
  using mapped_type = T;
  using value_type = std::pair<key_type, mapped_type>;
  using size_type = std::size_t;
  using const_iterator = const value_type*;

public:
  // Host methods
  bool finalized() const { return finalized_; }
  void insert(value_type value) { items_.push_back(value); }
  inline void erase(key_type key);
  inline void reserve(size_type count);
  inline void clear();
  inline void finalize();
  inline void copy_to_device();

  // Host/device methods
  bool empty() const { return items_.empty(); }
  size_type size() const { return items_.size(); }
  const_iterator begin() const { return items_.begin(); }
  const_iterator end() const { return items_.end(); }
  size_type bucket_count() const { return buckets_.size(); }
  inline size_type bucket_size(size_type idx) const;
  inline double load_factor() const;
  inline const T& operator[](key_type key) const;
  inline const_iterator find(key_type key) const;
  inline bool exists(key_type key) const;

private:
  // Data members
  vector<value_type> items_; //!< Key/value pairs
  vector<std::pair<size_type, size_type>>
    buckets_;              //!< Bucket start and end indices
  bool finalized_ {false}; //!< Whether hash table has been constructed
  Hash hash_;              //!< Hash function

  // Helper methods
  size_type bucket_index(const key_type& key) const;
  size_type calc_num_buckets(size_type count) const;

  // Sorting function for hashes
  struct HashComparator {
    const Hash& hash_;
    const int num_buckets_;

    HashComparator(const Hash& hash, int num_buckets)
      : hash_(hash), num_buckets_(num_buckets)
    {}

    bool operator()(const value_type& left, const value_type& right) const
    {
      auto lhash = hash_(left.first) % num_buckets_;
      auto rhash = hash_(right.first) % num_buckets_;

      if (lhash == rhash) {
        // Sort by key inside the hash
        return left.first < right.first;
      }
      return lhash < rhash;
    }
  };
};

//==========================================================================
// Inline definitions
//==========================================================================

template<typename Key, typename T, typename Hash>
auto static_map<Key, T, Hash>::bucket_size(size_type idx) const -> size_type
{
  const auto& bounds = buckets_[idx];
  return bounds.second - bounds.first;
}

template<typename Key, typename T, typename Hash>
double static_map<Key, T, Hash>::load_factor() const
{
  return static_cast<double>(this->size()) / this->bucket_count();
}

template<typename Key, typename T, typename Hash>
void static_map<Key, T, Hash>::erase(key_type key)
{
  // Search for the key
  auto iter = std::find_if(items_.begin(), items_.end(),
    [key](const value_type& kv) { return kv.first == key; });

  // Delete the item
  if (iter != items_.end()) {
    items_.erase(iter);
  }
}

template<typename Key, typename T, typename Hash>
void static_map<Key, T, Hash>::reserve(size_type count)
{
  items_.reserve(count);
  buckets_.reserve(this->calc_num_buckets(count));
}

template<typename Key, typename T, typename Hash>
void static_map<Key, T, Hash>::clear()
{
  items_.clear();
  buckets_.clear();
}

template<typename Key, typename T, typename Hash>
void static_map<Key, T, Hash>::finalize()
{
  // Empty map: nothing to do
  if (this->empty()) {
    finalized_ = true;
    return;
  }

  // Create empty buckets (begin == end)
  buckets_.resize(this->calc_num_buckets(this->size()), {0, 0});

  // Sort the key/value pairs according to hash, preserving the relative order
  // of duplicate items. This ensures that if elements with identical keys
  // exist, the first one inserted will be found. Note, however, that duplicates
  // are *not* removed from the vector -- it is up to the caller to avoid
  // inserting them unnecessarily.
  std::stable_sort(
    items_.begin(), items_.end(), HashComparator(hash_, this->bucket_count()));

  // Put the key/value pairs into consecutive buckets according to hash,
  // while checking for uniqueness by comparing against previous key
  auto prev_hash = this->bucket_index(items_.front().first);
  for (auto i = 1; i < this->size(); ++i) {
    const auto& key = items_[i].first;
    const auto hash = this->bucket_index(key);
    if (hash == prev_hash) {
      continue;
    }

    // If the hash is different, the last bucket has ended and the next
    // bucket has begun
    buckets_[prev_hash].second = i;
    buckets_[hash].first = i;
    prev_hash = hash;
  }
  // Add the final bucket 'end' index
  buckets_[prev_hash].second = this->size();
  finalized_ = true;
}

template<typename Key, typename T, typename Hash>
void static_map<Key, T, Hash>::copy_to_device()
{
    items_.copy_to_device();
    buckets_.copy_to_device();
}

template<typename Key, typename T, typename Hash>
const T& static_map<Key, T, Hash>::operator[](key_type key) const
{
  // Map must be finalized and key *must* exist in the map; for efficiency no
  // search is performed to see if the key has been inserted
  auto idx = this->bucket_index(key);
  const_iterator iter = items_.begin() + buckets_[idx].first;
  do {
    if (key == iter->first) {
      return iter->second;
    }
    ++iter;
  } while (true);
}

template<typename Key, typename T, typename Hash>
typename static_map<Key, T, Hash>::const_iterator
static_map<Key, T, Hash>::find(key_type key) const
{
  // Map must be finalized first!
  auto idx = this->bucket_index(key);
  const auto bounds = buckets_[idx];

  const auto last = items_.begin() + bounds.second;
  for (auto iter = items_.begin() + bounds.first; iter != last; ++iter) {
    if (key == iter->first) {
      return iter;
    }
  }
  // Not found
  return items_.end();
}

template<typename Key, typename T, typename Hash>
bool static_map<Key, T, Hash>::exists(key_type key) const
{
  // Map must be finalized first!
  return this->find(key) != items_.end();
}

template<typename Key, typename T, typename Hash>
auto static_map<Key, T, Hash>::bucket_index(const key_type& key) const
  -> size_type
{
  return hash_(key) % this->bucket_count();
}

} // namespace openmc

#endif
