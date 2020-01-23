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
// thread safe manner by use of atomics. 
template <typename T> 
class SharedArray { 

private: 
  //==========================================================================
  // Data members

  std::unique_ptr<T[]> data_; //!< A pointer to hold the data
  int64_t size_ {0}; //!< The current size of the shared array. 
  int64_t capacity_ {0}; //!< The maximum capacity of the shared array.

public: 
  //==========================================================================
  // Constructors

  // Creates an empty shared array
  SharedArray() = default;

  // Creates a shared array of desired capacity with zero size.
  SharedArray(int64_t capacity) : capacity_(capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
  }

  //==========================================================================
  // Methods and Accessors

  // Array accessor
  T& operator[](int64_t i) {return data_[i];}

  // Allocates space for the shared array to hold the indicated capacity
  void reserve(int64_t capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
    capacity_ = capacity;
  }

  // Increases the size of the SharedArray by one and returns an index to the
  // last element of the array.
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

  // Free's any space that was allocated to the shared array and resets
  // size and capacity to zero.
  void clear()
  {
    data_.reset();
    size_ = 0;
    capacity_ = 0;
  }

  // Size getter
  int64_t size() {return size_;}

  // Sets the size of the shared array. This is useful in cases where
  // we want to write to the shared array in a non-thread safe manner.
  int64_t resize(int64_t size) {size_ = size;}

  // Capacity getter
  int64_t capacity() {return capacity_;}

  // Returns a pointer to the data allocation
  T* data() {return data_.get();}

}; 

} // namespace openmc

#endif // OPENMC_SHARED_ARRAY_H
