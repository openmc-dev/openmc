#ifndef OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H
#define OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H

#include "openmc/openmp_interface.h"

#include <memory>

namespace openmc {

template<typename KeyType, typename ValueType, typename HashFunctor>
class ParallelMap {

  struct Bucket {
    OpenMPMutex lock_;
    std::unordered_map<KeyType, std::unique_ptr<ValueType>, HashFunctor> map_;
  };

public:
  ParallelMap(int n_buckets = 1000) : buckets_(n_buckets) {}

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

  ValueType* emplace(KeyType& key, ValueType& value)
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

private:
  Bucket& get_bucket(const KeyType& key)
  {
    return buckets_[hash(key) % buckets_.size()];
  }

  HashFunctor hash;
  vector<Bucket> buckets_;
};

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_PARALLEL_HASH_MAP_H