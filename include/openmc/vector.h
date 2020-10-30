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
#include <initializer_list>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace openmc {

/**
 * This is an implementation of a C++ standard vector,
 * with some extra cool features like thread safe append
 * on both CPU and GPU. It also compiles functions you may
 * hope to run on CUDA-capable devices into device code.
 *
 * Importantly, this uses CUDA unified memory. It's set up so
 * that all memory control should be done from the host, and
 * provides methods for data access that work on the device.
 * I'm not aware of any CUDA libraries which provide something
 * like this. It makes OpenMC GPU compatibility fit smoothly.
 *
 * By default, the vector size type is unsigned int as opposed
 * to the typical std::size_t, because we are _never_ going to
 * have vectors of lengths greater than a billion in Monte Carlo
 * apps. I would not template this on SizeType, but xtensor seems
 * to require size() to return std::size_t in order to properly
 * interface this class with it. As such, this class can be templated
 * to use std::size_t as its size type in the cases where that's
 * required.
 */
#ifdef __CUDACC__

template<typename T, typename Alloc = UnifiedAllocator<T>>
class vector {
public:
  using value_type = T;
  using size_type = unsigned int;
  using difference_type = int;
  using reference = T&;
  using const_reference = T const&;
  using pointer = typename Alloc::pointer;
  using iterator = T*;
  using const_iterator = T const*;
  using const_pointer = typename Alloc::const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using allocator_type = Alloc;

private:
  pointer begin_;
  unsigned int size_;
  unsigned int capacity_;
  Alloc alloc_;

public:
  __host__ void grow()
  {
    pointer old_begin = begin_;
    auto old_capacity = capacity_;
    capacity_ = capacity_ ? capacity_ * 2 : 1;
    begin_ = alloc_.allocate(capacity_);
    memcpy(begin_, old_begin, sizeof(T) * size_);
    if (old_capacity)
      alloc_.deallocate(old_begin, old_capacity);
  }

  __host__ vector() : begin_(nullptr), size_(0), capacity_(0) {}

  // Construct, default-initializing elements in a vector of
  // length new_size.
  __host__ vector(size_type new_size, const T& value = T())
  {
    begin_ = alloc_.allocate(new_size);
    capacity_ = new_size;
    for (size_type i = 0; i < new_size; ++i)
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
    for (size_type i = 0; i < size_; ++i)
      new (begin_ + i) T(copy_from[i]);
  }

  __host__ vector& operator=(vector&& move_from)
  {
    begin_ = std::move(move_from.begin_);
    size_ = std::move(move_from.size_);
    capacity_ = std::move(move_from.capacity_);
    alloc_ = std::move(move_from.alloc_);

    // Remove responsibity of the RHS for the memory it held
    move_from.begin_ = nullptr;
    move_from.size_ = 0;
    move_from.capacity_ = 0;

    return *this;
  }
  // Allow construction from vector
  __host__ vector& operator=(vector<T, Alloc> const& copy_from)
  {
    begin_ = alloc_.allocate(copy_from.size());
    size_ = copy_from.size();
    capacity_ = size_;
    alloc_ = copy_from.alloc_;
    for (size_type i = 0; i < size_; ++i)
      new (begin_ + i) T(copy_from[i]);
    return *this;
  }

  // Constructing from initializer lists
  __host__ vector(std::initializer_list<T> list)
    : begin_(nullptr), size_(0), capacity_(0)
  {
    resize(list.size());
    size_type i = 0;
    for (auto x : list)
      begin_[i++] = x;
  }

  // The move constructor may need to run on the device in the case of
  // construction of polymorphic objects living on GPU that contain vectors.
#pragma hd_warning_disable
  __host__ __device__ vector(vector&& move_from)
    : begin_(std::move(move_from.begin_)), size_(std::move(move_from.size_)),
      capacity_(std::move(move_from.capacity_)),
      alloc_(std::move(move_from.alloc_))
  {
    // Remove responsibility of the other one for the memory it held
    move_from.begin_ = nullptr;
    move_from.size_ = 0;
    move_from.capacity_ = 0;
  }

  __host__ ~vector()
  {
    for (size_type i = 0; i < size_; ++i)
      (begin_ + i)->~T();
    if (capacity_)
      alloc_.deallocate(begin_, capacity_);
  }

  // Construct from iterators
  __host__ vector(const_iterator begin, const_iterator end)
    : begin_(nullptr), size_(end - begin), capacity_(end - begin)
  {
    // If bad ordering, create vector of length zero
    if (end < begin) {
      size_ = 0;
      capacity_ = 0;
    }

    begin_ = alloc_.allocate(size_);
    for (size_type i = 0; i < size_; ++i)
      begin_[i] = begin[i];
  }

  // Construction from a a range of things that are convertible to value_type
  template<typename ConvertFrom,
    typename =
      std::enable_if_t<std::is_convertible<ConvertFrom, value_type>::value>>
  __host__ vector(const ConvertFrom* begin, const ConvertFrom* end)
    : begin_(nullptr), size_(end - begin), capacity_(end - begin)
  {
    // If bad ordering, create vector of length zero
    if (end < begin) {
      size_ = 0;
      capacity_ = 0;
    }
    begin_ = alloc_.allocate(size_);
    for (size_type i = 0; i < size_; ++i)
      new (begin_ + i) value_type(begin[i]);
  }

  // Conversion assignment
  template<typename ConvertFrom, typename Alloc2,
    typename =
      std::enable_if_t<std::is_convertible<ConvertFrom, value_type>::value>>
  __host__ vector& operator=(vector<ConvertFrom, Alloc2> const& rhs)
  {
    resize(rhs.size());
    for (size_type i = 0; i < rhs.size(); ++i)
      new (begin_ + i) value_type(rhs[i]);
    return *this;
  }

  __host__ void reserve(size_type new_size)
  {
    pointer old_begin = begin_;
    auto old_capacity = capacity_;
    capacity_ = new_size;
    begin_ = alloc_.allocate(capacity_);
    size_type num_to_copy = std::min(new_size, size_);
    memcpy(begin_, old_begin, num_to_copy * sizeof(T));
    if (old_capacity)
      alloc_.deallocate(old_begin, old_capacity);
  }

  __host__ void resize(size_type new_size)
  {
    if (new_size > capacity_)
      reserve(new_size);
    // Default initialize new things:
    for (size_type i = size_; i < new_size; ++i)
      begin_[i] = T();
    size_ = new_size;
  }
  __host__ void resize(size_type new_size, T const& default_value)
  {
    if (new_size > capacity_)
      reserve(new_size);
    // set new things to be default_value
    for (size_type i = size_; i < new_size; ++i)
      begin_[i] = default_value;
    size_ = new_size;
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
    // use copy constructor, NOT copy assignment which will deref the invalid
    // thing at begin_+size_
    new (begin_ + size_++) T(value);
  }

  __host__ size_type thread_safe_append(const T& value)
  {
    // This is what to do if acting on GPU
    size_type idx;
#ifdef __CUDA_ARCH__
    idx = atomicInc(&size_, capacity_);
#else
// else, do this on CPU
#pragma omp atomic capture
    idx = size_++;
    // For the applications here, it makes more sense to just do nothing
    if (idx >= capacity_) {
#pragma omp atomic write
      size_ = capacity_;
      idx = size_ - 1;
    }
#endif

    // Call copy constructor in the right place
    new (begin_ + idx) T(value);

    return idx;
  }

  __host__ void push_back(T&& value) { emplace_back(std::move(value)); }

  __host__ __device__ T& operator[](size_type n) { return *(begin_ + n); }
  __host__ __device__ T const& operator[](size_type n) const
  {
    return *(begin_ + n);
  }
  __host__ T& at(size_type n)
  {
    if (n >= size_) {
      throw std::out_of_range("openmc::vector::at() out of range!");
    }
    return *(begin_ + n);
  }
  __host__ T const& at(size_type n) const
  {
    if (n >= size_) {
      throw std::out_of_range("openmc::vector::at() out of range!");
    }
    return *(begin_ + n);
  }
  __host__ iterator insert(iterator pos, const T& value)
  {
    // Must handle edge case where vector has nothing in it yet
    size_type indx = pos - begin_;

    // Make space for the new entry if necessary
    if (size_ == capacity_)
      grow();

    for (size_type i = size_; i > indx; i--)
      begin_[i] = begin_[i - 1];
    begin_[indx] = value;
    size_++;
    return begin_ + indx;
  }

  // This one is pretty inefficient and could certainly be improved
  template<typename Iterator>
  __host__ void insert(iterator pos, Iterator first, Iterator last)
  {
    size_type indx = pos - begin_;
    // Have to recompute iterator since begin_ may be changing
    while (first != last)
      insert(begin_ + indx++, *(first++));
  }

  __host__ iterator erase(iterator pos)
  {
    size_type indx = pos - begin_;
    size_--;
    for (size_type i = indx; i < size_; ++i)
      begin_[i] = begin_[i + 1];
    return begin_ + indx;
  }

  __host__ __device__ bool empty() const { return size_ == 0; }
  __host__ __device__ pointer data() const { return begin_; }
  __host__ __device__ size_type size() const { return size_; }
  __host__ __device__ size_type capacity() const { return capacity_; }
  __host__ __device__ iterator begin() { return begin_; }
  __host__ __device__ iterator end() { return begin_ + size_; }
  __host__ __device__ iterator begin() const { return begin_; }
  __host__ __device__ iterator end() const { return begin_ + size_; }
  __host__ __device__ const_iterator cbegin() const { return begin_; }
  __host__ __device__ const_iterator cend() const { return begin_ + size_; }
  __host__ reverse_iterator rbegin() const
  {
    return std::make_reverse_iterator(begin_ + size_);
  }
  __host__ reverse_iterator rend() const
  {
    return std::make_reverse_iterator(begin_);
  }
  __host__ const_reverse_iterator crbegin() const
  {
    return std::make_reverse_iterator(begin_ + size_);
  }
  __host__ const_reverse_iterator crend() const
  {
    return std::make_reverse_iterator(begin_);
  }
  __host__ __device__ void pop_back() { size_--; }
  __host__ __device__ const_reference back() const { return begin_[size_ - 1]; }
  __host__ __device__ const_reference front() const { return begin_[0]; }
  __host__ __device__ reference back() { return begin_[size_ - 1]; }
  __host__ __device__ reference front() { return begin_[0]; }

  // Since in OpenMC, it is quite common to clear a vector and expect to
  // start pushing back to it again immediately after. So, it doesn't make
  // sense to deallocate memory here.
  __host__ void clear()
  {
    size_ = 0;
    for (size_type i = 0; i < size_; ++i)
      (begin_ + i)->~T();
  }

  // Enable host-device synchronization functions if replicated memory is being
  // used
  template<typename U = pointer>
  std::enable_if_t<std::is_same<U, DualPointer<T>>::value> syncToDevice()
  {
    cudaMemcpy(begin_.device_pointer(), begin_.host_pointer(),
      sizeof(value_type) * size_, cudaMemcpyHostToDevice);
  }
  template<typename U = pointer>
  std::enable_if_t<std::is_same<U, DualPointer<T>>::value> syncToHost()
  {
    cudaMemcpy(begin_.host_pointer(), begin_.device_pointer(),
      sizeof(value_type) * size_, cudaMemcpyDeviceToHost);
  }
};

template<typename T, typename Alloc>
__host__ __device__ bool operator==(
  const vector<T, Alloc>& first, const vector<T, Alloc>& second)
{
  if (first.size() != second.size())
    return false;
  for (typename vector<T, Alloc>::size_type i = 0; i < first.size(); ++i)
    if (first[i] != second[i])
      return false;
  return true;
}

#else

using std::vector;

#endif

} // namespace openmc
#endif // OPENMC_VECTOR_H
