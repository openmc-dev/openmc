#ifndef OPENMC_SHARED_ARRAY_H
#define OPENMC_SHARED_ARRAY_H

//! \file shared_array.h
//! \brief Shared array data structure

#include<memory>


namespace openmc {

//==============================================================================
// Class declarations
//==============================================================================

template <typename T> 
class SharedArray { 

private: 
  //==========================================================================
  // Data members
  std::unique_ptr<T[]> data_;
  int64_t size_ {0};
  int64_t capacity_ {0};

public: 
  //==========================================================================
  // Constructors
  
  SharedArray() = default;
  
  SharedArray(int64_t capacity) : capacity_(capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
  }
  
  //==========================================================================
  // Methods and Accessors
  
  T& operator[](int64_t i) {return data_[i];}

  void reserve(int64_t capacity)
  {
    data_ = std::make_unique<T[]>(capacity);
    capacity_ = capacity;
  }

  //! Increases the size of the SharedArray by one and returns an index to the
  //! last element of the array.
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

  void clear()
  {
    data_.reset();
    size_ = 0;
    capacity_ = 0;
  }

  int64_t size() {return size_;}
  
  int64_t resize(int64_t size) {size_ = size;}

  int64_t capacity() {return capacity_;}

  T* data() {return data_.get();}

}; 

} // namespace openmc

#endif // OPENMC_SHARED_ARRAY_H
