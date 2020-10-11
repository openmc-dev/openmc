#ifndef OPENMC_VECTOR_H
#define OPENMC_VECTOR_H

/*
 * In an implementation of OpenMC that offloads computations to an accelerator,
 * we may need to provide replacements for standard library containers and
 * algorithms that have no native implementations on the device of interest.
 * Because some developers are currently in the process of creating such code,
 * introducing the below typedef lessens the amount of rebase conflicts that
 * happen as they rebase their code on OpenMC's develop branch.
 *
 */

#include "openmc/memory.h"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace openmc {

/**
 * This is an implementation of a C++ standard vector.
 *
 * Importantly, this uses CUDA unified memory. It's set up so
 * that all memory control should be done from the host, and
 * provides methods for data access that work on the device.
 * I'm not aware of any CUDA libraries which provide something
 * like this. It makes OpenMC GPU compatibility fit smoothly.
 */
#ifdef __CUDACC__

template<class T, class Alloc = UnifiedAllocator<T>>
class vector {
private:
  T* begin_;
  std::size_t size_;
  std::size_t capacity_;
  Alloc alloc_;

  __host__ void grow()
  {
    T* old_begin = begin_;
    auto old_capacity = capacity_;
    // Multiply capacity by two if capacity is nonzero
    capacity_ = capacity_ ? capacity_ * 2 : 1;
    begin_ = alloc_.allocate(capacity_);
    memcpy(begin_, old_begin, size_ * sizeof(T));
    alloc_.deallocate(old_begin, old_capacity);
  }

public:
  __host__ vector() : begin_(nullptr), size_(0), capacity_(0) {}

  // Construct, default-initializing elements in a vector of
  // length new_size.
  __host__ vector(std::size_t new_size, const T& value = T())
  {
    begin_ = alloc_.allocate(new_size);
    capacity_ = new_size;
    for (std::size_t i = 0; i < new_size; ++i)
      new (begin_ + i) T(value);
    size_ = new_size;
  }

  // The copy constructor should create a totally independent deep copy
  __host__ vector(vector const& copy_from)
  {
    begin_ = alloc_.allocate(copy_from.size());
    size_ = copy_from.size();
    capacity_ = size_;
    alloc_ = copy_from.alloc_;
  }
  __host__ __device__ vector& operator=(vector&& copy_from)
  {
    begin_ = std::move(copy_from.begin_);
    size_ = std::move(copy_from.size_);
    capacity_ = std::move(copy_from.capacity_);
    alloc_ = std::move(copy_from.alloc_);
    return *this;
  }

  // The move constructor may need to run on the device in the case of
  // construction of polymorphic objects living on GPU that contain vectors.
  __host__ __device__ vector(vector&& move_from)
    : begin_(std::move(move_from.begin_)), size_(std::move(move_from.size_)),
      capacity_(std::move(move_from.capacity_)),
      alloc_(std::move(move_from.alloc_))
  {}

  __host__ ~vector()
  {
    for (std::size_t i = 0; i < size_; ++i)
      (begin_ + i)->~T();
    alloc_.deallocate(begin_, capacity_);
  }

  // Construct from iterators
  __host__ vector(const T* begin, const T* end)
    : begin_(nullptr), size_(end - begin), capacity_(end - begin)
  {
    // If bad ordering, create vector of length zero
    if (end < begin) {
      size_ = 0;
      capacity_ = 0;
    }

    begin_ = alloc_.allocate(size_);
    for (std::size_t i = 0; i < size_; ++i)
      begin_[i] = begin[i];
  }

  __host__ void reserve(std::size_t new_size)
  {
    T* old_begin = begin_;
    auto old_capacity = capacity_;
    capacity_ = new_size;
    begin_ = alloc_.allocate(capacity_);
    memcpy(begin_, old_begin, size_ * sizeof(T));
    alloc_.deallocate(old_begin, old_capacity);
  }

  __host__ void resize(std::size_t new_size)
  {
    reserve(new_size);
    // Default initialize new things:
    for (std::size_t i = size_; i < new_size; ++i)
      begin_[i] = T();
  }
  __host__ void resize(std::size_t new_size, T const& default_value)
  {
    reserve(new_size);
    // set new things to be default_value
    for (std::size_t i = size_; i < new_size; ++i)
      begin_[i] = default_value;
  }
  __host__ void shrink_to_fit() { reserve(size_); }

  template<class... Args>
  __host__ void emplace_back(Args&&... args)
  {
    if (size_ == capacity_)
      grow();
    new (begin_ + size_++) T(std::forward<Args>(args)...);
  }
  __host__ void push_back(const T& value)
  {
    if (size_ == capacity_)
      grow();
    begin_[size_++] = value;
  }
  __host__ void push_back(T&& value) { emplace_back(std::move(value)); }

  __host__ __device__ T& operator[](std::size_t n) { return *(begin_ + n); }
  __host__ __device__ T const& operator[](std::size_t n) const
  {
    return *(begin_ + n);
  }

  __host__ __device__ bool empty() const { return size_ == 0; }
  __host__ __device__ T* data() const { return begin_; }
  __host__ __device__ std::size_t size() const { return size_; }
  __host__ __device__ T* begin() { return begin_; }
  __host__ __device__ T* end() { return begin_ + size_; }

  __host__ void clear()
  {
    capacity_ = 0;
    size_ = 0;
    alloc_.deallocate(begin_, capacity_);
    begin_ = nullptr;
  }
};

#else

using std::vector;

#endif

} // namespace openmc
#endif // OPENMC_VECTOR_H
