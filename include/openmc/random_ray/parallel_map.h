#ifndef OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H
#define OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H

#include "openmc/openmp_interface.h"

#include <memory>

namespace openmc {

/*
 * The ParallelMap class allows for threadsafe access to a map-like data
 * structure. It is implemented as a hash table with a fixed number of buckets,
 * each of which contains a mutex lock and an unordered_map. The class provides
 * a subset of the functionality of std::unordered_map. Users must first lock
 * the object (using the key) before accessing or modifying the map. The object
 * is locked by bucket, allowing for multiple threads to manipulate different
 * keys simultaneously, though sometimes threads will need to wait if keys
 * happen to be in the same bucket.
 */

template<typename KeyType, typename ValueType, typename HashFunctor>
class ParallelMap {

  //----------------------------------------------------------------------------
  // Helper structes and classes

  struct Bucket {
    OpenMPMutex lock_;
    std::unordered_map<KeyType, std::unique_ptr<ValueType>, HashFunctor> map_;
  };

  // The iterator yields a pair: (const KeyType&, ValueType&)
  class iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<const KeyType&, ValueType&>;
    using difference_type = std::ptrdiff_t;
    using pointer = void; // Not providing pointer semantics.
    using reference = value_type;

    iterator(std::vector<Bucket>* buckets, std::size_t bucket_index,
      typename std::unordered_map<KeyType, std::unique_ptr<ValueType>,
        HashFunctor>::iterator inner_it)
      : buckets_(buckets), bucket_index_(bucket_index), inner_it_(inner_it)
    {
      // Advance to the first valid element if necessary.
      advance_to_valid();
    }

    // Dereference returns a pair of (key, value).
    reference operator*() const
    {
      return {inner_it_->first, *inner_it_->second};
    }

    iterator& operator++()
    {
      ++inner_it_;
      advance_to_valid();
      return *this;
    }

    iterator operator++(int)
    {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const iterator& other) const
    {
      // Two iterators are equal if they refer to the same bucket vector and are
      // both at end, or if they have the same bucket index and inner iterator.
      return buckets_ == other.buckets_ &&
             bucket_index_ == other.bucket_index_ &&
             (bucket_index_ == buckets_->size() ||
               inner_it_ == other.inner_it_);
    }

    bool operator!=(const iterator& other) const { return !(*this == other); }

  private:
    // Helper function: if we are at the end of the current bucket, advance to
    // the next non-empty bucket.
    void advance_to_valid()
    {
      while (bucket_index_ < buckets_->size() &&
             inner_it_ == (*buckets_)[bucket_index_].map_.end()) {
        ++bucket_index_;
        if (bucket_index_ < buckets_->size())
          inner_it_ = (*buckets_)[bucket_index_].map_.begin();
      }
    }

    std::vector<Bucket>* buckets_;
    std::size_t bucket_index_;
    typename std::unordered_map<KeyType, std::unique_ptr<ValueType>,
      HashFunctor>::iterator inner_it_;
  };

public:
  //----------------------------------------------------------------------------
  // Constructor
  ParallelMap(int n_buckets = 1000) : buckets_(n_buckets) {}

  //----------------------------------------------------------------------------
  // Public Methods
  void lock(const KeyType& key)
  {
    Bucket& bucket = get_bucket(key);
    bucket.lock_.lock();
  }

  void unlock(const KeyType& key)
  {
    Bucket& bucket = get_bucket(key);
    bucket.lock_.unlock();
  }

  void clear()
  {
    for (auto& bucket : buckets_) {
      bucket.map_.clear();
    }
  }

  bool contains(const KeyType& key)
  {
    Bucket& bucket = get_bucket(key);
    // C++20
    // return bucket.map_.contains(key);
    return bucket.map_.find(key) != bucket.map_.end();
  }

  ValueType& operator[](const KeyType& key)
  {
    Bucket& bucket = get_bucket(key);
    return *bucket.map_[key].get();
  }

  ValueType* emplace(KeyType& key, const ValueType& value)
  {
    Bucket& bucket = get_bucket(key);
    // Attempt to emplace the new element into the unordered_map within the
    auto result = bucket.map_.emplace(key, std::make_unique<ValueType>(value));
    auto it = result.first;
    return it->second.get();
  }

  // Copies everything into a vector, clears all buckets
  int64_t move_contents_into_vector(vector<ValueType>& v)
  {
    int64_t n_new = 0;
    for (auto& bucket : buckets_) {
      for (auto& pair : bucket.map_) {
        v.push_back(*pair.second.get());
        n_new++;
      }
      bucket.map_.clear();
    }
    return n_new;
  }

  // Return iterator to first element.
  iterator begin()
  {
    std::size_t bucket_index = 0;
    auto inner_it = buckets_.empty()
                      ? typename std::unordered_map<KeyType,
                          std::unique_ptr<ValueType>, HashFunctor>::iterator()
                      : buckets_[0].map_.begin();
    return iterator(&buckets_, bucket_index, inner_it);
  }

  // Return iterator to one-past-last element.
  iterator end()
  {
    // End is signaled by bucket_index_ equal to buckets_.size()
    return iterator(&buckets_, buckets_.size(),
      typename std::unordered_map<KeyType, std::unique_ptr<ValueType>,
        HashFunctor>::iterator());
  }

private:
  //----------------------------------------------------------------------------
  // Private Methods
  Bucket& get_bucket(const KeyType& key)
  {
    return buckets_[hash(key) % buckets_.size()];
  }

  //----------------------------------------------------------------------------
  // Private Data Fields
  HashFunctor hash;
  vector<Bucket> buckets_;
};

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H