#ifndef OPENMC_SHARED_ARRAY_H
#define OPENMC_SHARED_ARRAY_H

//! \file shared_array.h
//! \brief Shared array data structure

#include <memory>


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
template <typename T> 
class SharedArray { 

public: 
  //==========================================================================
  // Constructors

  //! Default constructor.
  SharedArray() = default;

  //! Construct a zero size container with space to hold capacity number of
  //! elements.
  //
  //! \param capacity The number of elements for the container to allocate
  //! space for
  SharedArray(int64_t capacity) : capacity_(capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
  }

  //==========================================================================
  // Methods and Accessors

  //! Return a reference to the element at specified location i. No bounds
  //! checking is performed.
  T& operator[](int64_t i) {return data_[i];}
  const T& operator[](int64_t i) const { return data_[i]; }

  //! Allocate space in the container for the specified number of elements.
  //! reserve() does not change the size of the container.
  //
  //! \param capacity The number of elements to allocate in the container
  void reserve(int64_t capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
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
    #pragma omp atomic capture
    idx = size_++;

    // Check that we haven't written off the end of the array
    if (idx >= capacity_) {
      #pragma omp atomic write
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

  //! Return the number of elements in the container
  int64_t size() {return size_;}

  //! Resize the container to contain a specified number of elements. This is
  //! useful in cases where the container is written to in a non-thread safe manner,
  //! where the internal size of the array needs to be manually updated.
  //
  //! \param size The new size of the container
  void resize(int64_t size) {size_ = size;}

  //! Return the number of elements that the container has currently allocated
  //! space for.
  int64_t capacity() {return capacity_;}

  //! Return pointer to the underlying array serving as element storage.
  T* data() {return data_.get();}
  const T* data() const {return data_.get();}

private: 
  //==========================================================================
  // Data members

  std::unique_ptr<T[]> data_; //!< An RAII handle to the elements
  int64_t size_ {0}; //!< The current number of elements 
  int64_t capacity_ {0}; //!< The total space allocated for elements

}; 

} // namespace openmc

#endif // OPENMC_SHARED_ARRAY_H
