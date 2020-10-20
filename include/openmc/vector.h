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

  // For stuff like std::string, you explicitly NEED to call
  // a copy assignment constructor. memcpy of std::string will fail.
  // template<typename U = T>
  // __host__ std::enable_if_t<std::is_copy_assignable<U>::value> grow()
  // {
  //   T* old_begin = begin_;
  //   auto old_capacity = capacity_;
  //   capacity_ = capacity_ ? capacity_ * 2 : 1;
  //   begin_ = alloc_.allocate(capacity_);
  //   // std::copy(old_begin, old_begin+size_, begin_);
  //   for (std::size_t i=0; i<size_; ++i)
  //     new(begin_+i) T(old_begin[i]);
  //   alloc_.deallocate(old_begin, old_capacity);
  // }
  __host__ void grow()
  {
    T* old_begin = begin_;
    auto old_capacity = capacity_;
    capacity_ = capacity_ ? capacity_ * 2 : 1;
    begin_ = alloc_.allocate(capacity_);
    memcpy(begin_, old_begin, sizeof(T) * size_);
    if (old_capacity)
      alloc_.deallocate(old_begin, old_capacity);
  }

  // // For stuff like unique_ptr _without_ copy assignment, just directly copy
  // template<typename U = T>
  // __host__ std::enable_if_t<!std::is_copy_assignable<U>::value> grow()
  // {
  //   T* old_begin = begin_;
  //   auto old_capacity = capacity_;
  //   capacity_ = capacity_ ? capacity_ * 2 : 1;
  //   begin_ = alloc_.allocate(capacity_);
  //   memcpy(begin_, old_begin, sizeof(T) * size_);
  //   alloc_.deallocate(old_begin, old_capacity);
  // }

public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::size_t;
  using reference = T&;
  using const_reference = T const&;
  using pointer = T*;
  using iterator = T*;
  using const_iterator = T const*;
  using const_pointer = T const*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using allocator_type = Alloc;

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
    for (std::size_t i = 0; i < size_; ++i)
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
  __host__ vector& operator=(vector<T> const& copy_from)
  {
    begin_ = alloc_.allocate(copy_from.size());
    size_ = copy_from.size();
    capacity_ = size_;
    alloc_ = copy_from.alloc_;
    for (std::size_t i = 0; i < size_; ++i)
      new (begin_ + i) T(copy_from[i]);
    return *this;
  }
  // template<typename = std::enable_if_t<
  // __host__ vector& operator=(std::vector<T>&& copy_from)
  // {
  //   // TODO this could be improved to steal the memory the other vector
  //   // is holding. However, this memory would not  be compatible with the
  //   // CUDA managed memory system unless that vector also used a
  //   // UnifiedAllocator, so this could be specialized for that special
  //   // case where the memory is known to be stealable. I am unsure how to
  //   // make the old vector not call delete on its memory when it's done
  //   // though.
  //   resize(copy_from.size());
  //   for (std::size_t i = 0; i < size_; ++i)
  //     begin_[i] = std::move(copy_from[i]);
  //   return *this;
  // }

  // Constructing from initializer lists
  __host__ vector(std::initializer_list<T> list)
    : begin_(nullptr), size_(0), capacity_(0)
  {
    resize(list.size());
    std::size_t i = 0;
    for (auto x : list)
      begin_[i++] = x;
  }

  // The move constructor may need to run on the device in the case of
  // construction of polymorphic objects living on GPU that contain vectors.
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
    for (std::size_t i = 0; i < size_; ++i)
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
    for (std::size_t i = 0; i < size_; ++i)
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
    for (std::size_t i = 0; i < size_; ++i)
      begin_[i] = begin[i];
  }

  // Conversion assignment
  template<typename ConvertFrom, typename Alloc2,
    typename =
      std::enable_if_t<std::is_convertible<ConvertFrom, value_type>::value>>
  __host__ vector& operator=(vector<ConvertFrom, Alloc2> const& rhs)
  {
    resize(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i)
      new (begin_ + i) value_type(rhs[i]);
    return *this;
  }

  __host__ void reserve(std::size_t new_size)
  {
    T* old_begin = begin_;
    auto old_capacity = capacity_;
    capacity_ = new_size;
    begin_ = alloc_.allocate(capacity_);
    std::size_t num_to_copy = std::min(new_size, size_);
    memcpy(begin_, old_begin, num_to_copy * sizeof(T));
    if (old_capacity)
      alloc_.deallocate(old_begin, old_capacity);
  }

  __host__ void resize(std::size_t new_size)
  {
    reserve(new_size);
    size_ = new_size;
    // Default initialize new things:
    for (std::size_t i = size_; i < new_size; ++i)
      begin_[i] = T();
  }
  __host__ void resize(std::size_t new_size, T const& default_value)
  {
    reserve(new_size);
    size_ = new_size;
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
    // use copy constructor, NOT copy assignment which will deref the invalid
    // thing at begin_+size_
    new (begin_ + size_++) T(value);
  }
  __host__ void push_back(T&& value) { emplace_back(std::move(value)); }

  __host__ __device__ T& operator[](std::size_t n) { return *(begin_ + n); }
  __host__ __device__ T const& operator[](std::size_t n) const
  {
    return *(begin_ + n);
  }
  __host__ T& at(std::size_t n)
  {
    if (n >= size_) {
      throw std::out_of_range("openmc::vector::at() out of range!");
    }
    return *(begin_ + n);
  }
  __host__ T const& at(std::size_t n) const
  {
    if (n >= size_) {
      throw std::out_of_range("openmc::vector::at() out of range!");
    }
    return *(begin_ + n);
  }
  __host__ iterator insert(iterator pos, const T& value)
  {
    // Make space for the new entry, and grow if necessary
    if (size_ == capacity_)
      grow();
    std::size_t indx = pos - begin_;
    for (std::size_t i = size_; i > indx; i--)
      begin_[i] = begin_[i - 1];
    begin_[indx] = value;
    size_++;
    return begin_ + indx;
  }

  template<typename Iterator>
  __host__ void insert(iterator pos, Iterator first, Iterator last)
  {
    while (first != last)
      insert(pos, *(first++));
  }

  __host__ iterator erase(iterator pos)
  {
    std::size_t indx = pos - begin_;
    size_--;
    for (std::size_t i = indx; i < size_; ++i)
      begin_[i] = begin_[i + 1];
    return begin_ + indx;
  }

  __host__ __device__ bool empty() const { return size_ == 0; }
  __host__ __device__ T* data() const { return begin_; }
  __host__ __device__ std::size_t size() const { return size_; }
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
  __host__ __device__ reference back() const { return begin_[size_ - 1]; }
  __host__ __device__ reference front() const { return begin_[0]; }

  __host__ void clear()
  {
    capacity_ = 0;
    size_ = 0;
    for (std::size_t i = 0; i < size_; ++i)
      (begin_ + i)->~T();
    alloc_.deallocate(begin_, capacity_);
    begin_ = nullptr;
  }
};

template<typename T, typename Alloc>
__host__ __device__ bool operator==(
  const vector<T, Alloc>& first, const vector<T, Alloc>& second)
{
  if (first.size() != second.size())
    return false;
  for (std::size_t i = 0; i < first.size(); ++i)
    if (first[i] != second[i])
      return false;
  return true;
}

#else

using vector;

#endif

} // namespace openmc
#endif // OPENMC_VECTOR_H
