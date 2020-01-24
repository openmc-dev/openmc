#ifndef OPENMC_SHARED_ARRAY_H
#define OPENMC_SHARED_ARRAY_H

//! \file shared_array.h
//! \brief Shared array data structure

#include<memory>


namespace openmc {

//==============================================================================
// Class declarations
//==============================================================================

// The SharedArray is an array that is capable of being appended to in an
// thread safe manner by use of atomics. It only provides protection for the
// use cases currently present in OpenMC. Namely, it covers the scenario where
// multiple threads are appending to an array, but no threads are reading from
// or operating on it in any other way at the same time. Multiple threads can
// call the thread_safe_append() function concurrently and store data to the
// object at the index returned from thread_safe_append() safely, but no other
// operations are protected.
template <typename T> 
class SharedArray { 

private: 
  //==========================================================================
  // Data members

  std::unique_ptr<T[]> data_; //!< A pointer to hold the data
  int64_t size_ {0}; //!< The current size of the SharedArray. 
  int64_t capacity_ {0}; //!< The maximum capacity of the SharedArray.

public: 
  //==========================================================================
  // Constructors

  //! Creates an empty SharedArray
  SharedArray() = default;

  //! Creates a SharedArray of desired capacity with zero size.
  //! \param capacity The desired capacity to allocate for the array
  SharedArray(int64_t capacity) : capacity_(capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
  }

  //==========================================================================
  // Methods and Accessors

  //! Array accessor
  T& operator[](int64_t i) {return data_[i];}

  //! Allocates space for the SharedArray
  //! \param capacity The number of elements to allocate in the array.
  void reserve(int64_t capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
    capacity_ = capacity;
  }

  //! Increases the size of the SharedArray by one and returns an index to the
  //! last element of the array. Also tests to enforce that the append
  //! operation does not read off the end of the array. In the event that this
  //! does happen, the size is set to be equal to the capacity, and -1 is 
  //! returned.
  //! \return The last index in the array, which is safe to write to. In the
  //! event that this index would be greater than what was allocated for the
  //! SharedArray, -1 is returned.
  int64_t thread_safe_append()
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

    return idx;
  }

  //! Frees any space that was allocated to the SharedArray and resets
  //! size and capacity to zero.
  void clear()
  {
    data_.reset();
    size_ = 0;
    capacity_ = 0;
  }

  //! Size getter function
  //! \return The current size of the SharedArray
  int64_t size() {return size_;}

  //! Sets the size of the SharedArray. This is useful in cases where
  //! we want to write to the SharedArray in a non-thread safe manner
  //! and need to update the internal size of the array after doing so.
  //! \param size The new size for the array
  void resize(int64_t size) {size_ = size;}

  //! Capacity getter functon
  //! \return The maximum allocated capacity for the SharedArray
  int64_t capacity() {return capacity_;}

  //! This function exposes the pointer that the SharedArray is protecting.
  //! \return The pointer to the data allocation
  T* data() {return data_.get();}

}; 

} // namespace openmc

#endif // OPENMC_SHARED_ARRAY_H
