#ifndef OPENMC_SHARED_ARRAY_H
#define OPENMC_SHARED_ARRAY_H

//! \file shared_array.h
//! \brief Shared array data structure

#include <algorithm> // for copy_n

#include "openmc/memory.h"

namespace openmc {

//==============================================================================
// Class declarations
//==============================================================================

// This container is an array that is capable of being appended to in an
// thread safe manner by use of atomics. It only provides protection for the
// use cases currently present in OpenMC. Namely, it covers the scenario where
// multiple threads are appending to an array, but no threads are reading from
// or operating on it in any other way at the same time. Multiple threads can
// call the thread_safe_append() function concurrently and store data to the
// object at the index returned from thread_safe_append() safely, but no other
// operations are protected.
template<typename T>
class SharedArray {

public:
  //==========================================================================
  // Constructors

  //! Default constructor.
  SharedArray() = default;

  //! Construct a container with `size` elements and capacity equal to `size`.
  //
  //! \param size The number of elements to allocate and initialize
  SharedArray(int64_t size) : size_(size), capacity_(size)
  {
    data_ = make_unique<T[]>(size);
  }

  //==========================================================================
  // Methods and Accessors

  //! Return a reference to the element at specified location i. No bounds
  //! checking is performed.
  T& operator[](int64_t i) { return data_[i]; }
  const T& operator[](int64_t i) const { return data_[i]; }

  //! Allocate space in the container for the specified number of elements.
  //! reserve() does not change the size of the container.
  //
  //! \param capacity The number of elements to allocate in the container
  void reserve(int64_t capacity)
  {
    data_ = make_unique<T[]>(capacity);
    capacity_ = capacity;
  }

  //! Increase the size of the container by one and append value to the
  //! array. Returns an index to the element of the array written to. Also
  //! tests to enforce that the append operation does not read off the end
  //! of the array. In the event that this does happen, set the size to be
  //! equal to the capacity and return -1.
  //
  //! \value The value of the element to append
  //! \return The index in the array written to. In the event that this
  //! index would be greater than what was allocated for the container,
  //! return -1.
  int64_t thread_safe_append(const T& value)
  {
    // Atomically capture the index we want to write to
    int64_t idx;
#pragma omp atomic capture seq_cst
    idx = size_++;

    // Check that we haven't written off the end of the array
    if (idx >= capacity_) {
#pragma omp atomic write seq_cst
      size_ = capacity_;
      return -1;
    }

    // Copy element value to the array
    data_[idx] = value;

    return idx;
  }

  //! Free any space that was allocated for the container. Set the
  //! container's size and capacity to 0.
  void clear()
  {
    data_.reset();
    size_ = 0;
    capacity_ = 0;
  }

  //! Push back an element to the array, with capacity and reallocation behavior
  //! as if this were a vector. This does not perform any thread safety checks.
  //! If the size exceeds the capacity, then the capacity will double just as
  //! with a vector. Data will be reallocated and moved to a new pointer and
  //! copied in before the new item is appended. Old data will be freed.
  void thread_unsafe_append(const T& value)
  {
    if (size_ == capacity_) {
      int64_t new_capacity = capacity_ == 0 ? 8 : 2 * capacity_;
      unique_ptr<T[]> new_data = make_unique<T[]>(new_capacity);
      std::copy_n(data_.get(), size_, new_data.get());
      data_ = std::move(new_data);
      capacity_ = new_capacity;
    }
    data_[size_++] = value;
  }

  //! Increase the size of the container by count elements without assigning
  //! values to the new elements. Existing elements are preserved if the
  //! container needs to grow. This does not perform any thread safety checks.
  //
  //! \param count The number of elements to append
  //! \return The starting index of the appended range
  int64_t extend_uninitialized(int64_t count)
  {
    int64_t offset = size_;
    int64_t new_size = size_ + count;
    if (new_size > capacity_) {
      int64_t new_capacity = capacity_ == 0 ? 8 : capacity_;
      while (new_capacity < new_size) {
        new_capacity *= 2;
      }
      unique_ptr<T[]> new_data = make_unique<T[]>(new_capacity);
      if (size_ > 0) {
        std::copy_n(data_.get(), size_, new_data.get());
      }
      data_ = std::move(new_data);
      capacity_ = new_capacity;
    }
    size_ = new_size;
    return offset;
  }

  //! Return the number of elements in the container
  int64_t size() { return size_; }
  int64_t size() const { return size_; }

  //! Resize the container to contain a specified number of elements. This is
  //! useful in cases where the container is written to in a non-thread safe
  //! manner, where the internal size of the array needs to be manually updated.
  //
  //! \param size The new size of the container
  void resize(int64_t size) { size_ = size; }

  //! Return whether the array is full
  bool full() const { return size_ == capacity_; }

  //! Return the number of elements that the container has currently allocated
  //! space for.
  int64_t capacity() { return capacity_; }

  //! Return pointer to the underlying array serving as element storage.
  T* data() { return data_.get(); }
  const T* data() const { return data_.get(); }

  //! Classic iterators
  T* begin() { return data_.get(); }
  const T* cbegin() const { return data_.get(); }
  T* end() { return data_.get() + size_; }
  const T* cend() const { return data_.get() + size_; }

private:
  //==========================================================================
  // Data members

  unique_ptr<T[]> data_; //!< An RAII handle to the elements
  int64_t size_ {0};     //!< The current number of elements
  int64_t capacity_ {0}; //!< The total space allocated for elements
};

} // namespace openmc

#endif // OPENMC_SHARED_ARRAY_H
