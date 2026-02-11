//! \file tensor.h
//! \brief Multi-dimensional tensor types for OpenMC.
//!
//! Provides Tensor<T> (dynamic-rank), Fixed2D<T,R,C> (stack-allocated),
//! and lightweight view types (View1D, ViewFlat).

#ifndef OPENMC_TENSOR_H
#define OPENMC_TENSOR_H

#include "openmc/vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <type_traits>

namespace openmc {
namespace tensor {

//==============================================================================
// Forward declarations
//==============================================================================

template<typename T>
class Tensor;

template<typename T, size_t R, size_t C>
class Fixed2D;

//==============================================================================
// Storage type mapping
//
// std::vector<bool> is a bit-packed specialization that returns proxy objects
// instead of real references, which breaks generic code.  storage_type_map
// redirects bool to unsigned char so that Tensor<bool> stores one byte per
// element with normal reference semantics.
//==============================================================================

template<typename T>
struct storage_type_map {
  using type = T;
};
template<>
struct storage_type_map<bool> {
  using type = unsigned char;
};
template<typename T>
using storage_type = typename storage_type_map<T>::type;

//==============================================================================
// View1D<T>: a read/write view of one row, column, or slice of a tensor.
//
// Holds a pointer, element count, and stride into the parent tensor's
// storage — no allocation or copy.
//==============================================================================

template<typename T>
class View1D {
public:
  using value_type = std::remove_const_t<T>;

  View1D(T* data, size_t size, size_t stride = 1)
    : data_(data), size_(size), stride_(stride)
  {}

  T& operator()(size_t i) { return data_[i * stride_]; }
  const T& operator()(size_t i) const { return data_[i * stride_]; }
  T& operator[](size_t i) { return data_[i * stride_]; }
  const T& operator[](size_t i) const { return data_[i * stride_]; }

  size_t size() const { return size_; }
  T* data() { return data_; }
  const T* data() const { return data_; }
  size_t stride() const { return stride_; }

  View1D<T> slice(size_t start, size_t end)
  {
    return {data_ + start * stride_, end - start, stride_};
  }
  View1D<const T> slice(size_t start, size_t end) const
  {
    return {data_ + start * stride_, end - start, stride_};
  }
  View1D<T> slice(size_t start)
  {
    return {data_ + start * stride_, size_ - start, stride_};
  }

  // Assignment from scalar
  template<typename U>
  auto operator=(U val) ->
    std::enable_if_t<std::is_arithmetic<U>::value, View1D&>
  {
    for (size_t i = 0; i < size_; ++i)
      data_[i * stride_] = val;
    return *this;
  }

  // Assignment from initializer_list
  View1D& operator=(std::initializer_list<value_type> vals)
  {
    auto it = vals.begin();
    for (size_t i = 0; i < size_ && it != vals.end(); ++i, ++it)
      data_[i * stride_] = *it;
    return *this;
  }

  // Assignment from Tensor (deferred, defined after Tensor)
  template<typename U>
  View1D& operator=(const Tensor<U>& other);

  // Compound assignment from Tensor (deferred)
  template<typename U>
  View1D& operator+=(const Tensor<U>& o);

  View1D& operator*=(value_type val)
  {
    for (size_t i = 0; i < size_; ++i)
      data_[i * stride_] *= val;
    return *this;
  }

  // Iterators
  class const_iterator {
    const T* ptr_;
    size_t stride_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_const_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    const_iterator(const T* ptr, size_t stride)
      : ptr_(ptr), stride_(stride)
    {}
    const T& operator*() const { return *ptr_; }
    const_iterator& operator++()
    {
      ptr_ += stride_;
      return *this;
    }
    const_iterator operator++(int)
    {
      auto tmp = *this;
      ptr_ += stride_;
      return tmp;
    }
    const_iterator& operator--()
    {
      ptr_ -= stride_;
      return *this;
    }
    const_iterator operator+(difference_type n) const
    {
      return const_iterator(ptr_ + n * stride_, stride_);
    }
    const_iterator operator-(difference_type n) const
    {
      return const_iterator(ptr_ - n * stride_, stride_);
    }
    difference_type operator-(const const_iterator& other) const
    {
      return (ptr_ - other.ptr_) / static_cast<difference_type>(stride_);
    }
    bool operator==(const const_iterator& other) const
    {
      return ptr_ == other.ptr_;
    }
    bool operator!=(const const_iterator& other) const
    {
      return ptr_ != other.ptr_;
    }
    bool operator<(const const_iterator& other) const
    {
      return ptr_ < other.ptr_;
    }
    bool operator>(const const_iterator& other) const
    {
      return ptr_ > other.ptr_;
    }
    bool operator<=(const const_iterator& other) const
    {
      return ptr_ <= other.ptr_;
    }
    bool operator>=(const const_iterator& other) const
    {
      return ptr_ >= other.ptr_;
    }
    const T& operator[](difference_type n) const
    {
      return *(ptr_ + n * stride_);
    }
    const_iterator& operator+=(difference_type n)
    {
      ptr_ += n * stride_;
      return *this;
    }
    const_iterator& operator-=(difference_type n)
    {
      ptr_ -= n * stride_;
      return *this;
    }
    friend const_iterator operator+(
      difference_type n, const const_iterator& it)
    {
      return it + n;
    }
  };

  class iterator {
    T* ptr_;
    size_t stride_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_const_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    iterator(T* ptr, size_t stride) : ptr_(ptr), stride_(stride) {}
    T& operator*() { return *ptr_; }
    iterator& operator++()
    {
      ptr_ += stride_;
      return *this;
    }
    iterator operator++(int)
    {
      auto tmp = *this;
      ptr_ += stride_;
      return tmp;
    }
    iterator& operator--()
    {
      ptr_ -= stride_;
      return *this;
    }
    iterator operator+(difference_type n) const
    {
      return iterator(ptr_ + n * stride_, stride_);
    }
    iterator operator-(difference_type n) const
    {
      return iterator(ptr_ - n * stride_, stride_);
    }
    difference_type operator-(const iterator& other) const
    {
      return (ptr_ - other.ptr_) / static_cast<difference_type>(stride_);
    }
    bool operator==(const iterator& other) const
    {
      return ptr_ == other.ptr_;
    }
    bool operator!=(const iterator& other) const
    {
      return ptr_ != other.ptr_;
    }
    bool operator<(const iterator& other) const
    {
      return ptr_ < other.ptr_;
    }
    T& operator[](difference_type n) { return *(ptr_ + n * stride_); }
    iterator& operator+=(difference_type n)
    {
      ptr_ += n * stride_;
      return *this;
    }
    friend iterator operator+(difference_type n, const iterator& it)
    {
      return it + n;
    }
  };

  iterator begin() { return iterator(data_, stride_); }
  iterator end() { return iterator(data_ + size_ * stride_, stride_); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }
  const_iterator cbegin() const { return const_iterator(data_, stride_); }
  const_iterator cend() const
  {
    return const_iterator(data_ + size_ * stride_, stride_);
  }

private:
  T* data_;
  size_t size_;
  size_t stride_;
};

//==============================================================================
// ViewFlat<T>: a flat view of all elements of a tensor.
//==============================================================================

template<typename T>
class ViewFlat {
public:
  ViewFlat(T* data, size_t size) : data_(data), size_(size) {}

  T& operator()(size_t i) { return data_[i]; }
  const T& operator()(size_t i) const { return data_[i]; }

  template<typename U>
  auto operator=(U val) ->
    std::enable_if_t<std::is_arithmetic<U>::value, ViewFlat&>
  {
    std::fill(data_, data_ + size_, static_cast<T>(val));
    return *this;
  }

  T* data() { return data_; }
  const T* data() const { return data_; }
  size_t size() const { return size_; }

  T* begin() { return data_; }
  T* end() { return data_ + size_; }
  const T* begin() const { return data_; }
  const T* end() const { return data_ + size_; }

private:
  T* data_;
  size_t size_;
};

//==============================================================================
// Tensor<T>: dynamic-rank N-dimensional tensor.
//
// Stores elements in a contiguous row-major vector<storage_type<T>>
// with a dynamic shape.
//==============================================================================

template<typename T>
class Tensor {
public:
  using value_type = T;
  using stored_type = storage_type<T>;
  using iterator = typename vector<stored_type>::iterator;
  using const_iterator = typename vector<stored_type>::const_iterator;

  //--------------------------------------------------------------------------
  // Constructors

  Tensor() = default;

  //! Construct with shape (uninitialized for arithmetic types via vector resize)
  explicit Tensor(vector<size_t> shape)
    : shape_(std::move(shape)), data_(compute_size())
  {}

  //! Construct with shape and fill value
  Tensor(vector<size_t> shape, T fill)
    : shape_(std::move(shape)), data_(compute_size(), fill)
  {}

  //! Construct from initializer_list shape
  explicit Tensor(std::initializer_list<size_t> shape)
    : shape_(shape), data_(compute_size())
  {}

  //! Construct from initializer_list shape with fill
  Tensor(std::initializer_list<size_t> shape, T fill)
    : shape_(shape), data_(compute_size(), fill)
  {}

  //! 1D copy from vector (disabled when T=size_t to avoid ambiguity)
  template<typename Dummy = T,
    typename = std::enable_if_t<!std::is_same<Dummy, size_t>::value>>
  explicit Tensor(const vector<T>& vec)
    : shape_({vec.size()}), data_(vec.begin(), vec.end())
  {}

  //! 1D copy from std::vector (disabled when T=size_t)
  template<typename Alloc, typename Dummy = T,
    typename = std::enable_if_t<!std::is_same<Dummy, size_t>::value>>
  explicit Tensor(const std::vector<T, Alloc>& vec)
    : shape_({vec.size()}), data_(vec.begin(), vec.end())
  {}

  //! 1D copy from raw pointer + count
  Tensor(const T* ptr, size_t count)
    : shape_({count}), data_(ptr, ptr + count)
  {}

  //! Copy from vector with explicit shape
  template<typename Alloc>
  Tensor(const std::vector<T, Alloc>& vec, vector<size_t> shape)
    : shape_(std::move(shape)), data_(vec.begin(), vec.end())
  {}

  //! Copy from View1D (makes a 1D tensor)
  template<typename U>
  Tensor(const View1D<U>& v)
    : shape_({v.size()})
  {
    data_.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i)
      data_[i] = v(i);
  }

  //! Cross-type copy constructor
  template<typename U,
    typename = std::enable_if_t<!std::is_same<U, T>::value>>
  Tensor(const Tensor<U>& other)
    : shape_(other.shape())
  {
    data_.resize(other.size());
    for (size_t i = 0; i < other.size(); ++i)
      data_[i] = static_cast<stored_type>(other.data()[i]);
  }

  //--------------------------------------------------------------------------
  // Assignment

  //! Cross-type assignment
  template<typename U,
    typename = std::enable_if_t<!std::is_same<U, T>::value>>
  Tensor& operator=(const Tensor<U>& other)
  {
    shape_ = other.shape();
    data_.resize(other.size());
    for (size_t i = 0; i < other.size(); ++i)
      data_[i] = static_cast<stored_type>(other.data()[i]);
    return *this;
  }

  //! Assignment from View1D
  Tensor& operator=(const View1D<T>& v)
  {
    shape_ = {v.size()};
    data_.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i)
      data_[i] = v(i);
    return *this;
  }

  //! Assignment from initializer_list of values (1D)
  template<typename Dummy = T,
    typename = std::enable_if_t<!std::is_same<Dummy, size_t>::value>>
  Tensor& operator=(std::initializer_list<T> vals)
  {
    shape_ = {vals.size()};
    data_.assign(vals.begin(), vals.end());
    return *this;
  }

  //--------------------------------------------------------------------------
  // Accessors

  stored_type* data() { return data_.data(); }
  const stored_type* data() const { return data_.data(); }
  size_t size() const { return data_.size(); }
  const vector<size_t>& shape() const { return shape_; }
  size_t shape(size_t dim) const {
    return dim < shape_.size() ? shape_[dim] : 0;
  }
  size_t ndim() const { return shape_.size(); }
  bool empty() const { return data_.empty(); }

  //--------------------------------------------------------------------------
  // Indexing (row-major)

  template<typename... Indices>
  stored_type& operator()(Indices... indices)
  {
    const size_t idx[] = {static_cast<size_t>(indices)...};
    size_t off = 0;
    for (size_t d = 0; d < sizeof...(Indices); ++d)
      off = off * shape_[d] + idx[d];
    return data_[off];
  }

  template<typename... Indices>
  const stored_type& operator()(Indices... indices) const
  {
    const size_t idx[] = {static_cast<size_t>(indices)...};
    size_t off = 0;
    for (size_t d = 0; d < sizeof...(Indices); ++d)
      off = off * shape_[d] + idx[d];
    return data_[off];
  }

  stored_type& operator[](size_t i) { return data_[i]; }
  const stored_type& operator[](size_t i) const { return data_[i]; }

  //--------------------------------------------------------------------------
  // Iterators

  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }
  const_iterator cbegin() const { return data_.cbegin(); }
  const_iterator cend() const { return data_.cend(); }

  //--------------------------------------------------------------------------
  // Mutation

  void resize(const vector<size_t>& shape)
  {
    shape_ = shape;
    data_.resize(compute_size());
  }

  void resize(std::initializer_list<size_t> shape)
  {
    shape_.assign(shape.begin(), shape.end());
    data_.resize(compute_size());
  }

  void resize(const vector<unsigned long long>& shape)
  {
    shape_.clear();
    for (auto d : shape)
      shape_.push_back(static_cast<size_t>(d));
    data_.resize(compute_size());
  }

  template<typename ShapeType>
  void reshape(const ShapeType& new_shape)
  {
    shape_.clear();
    for (auto d : new_shape)
      shape_.push_back(static_cast<size_t>(d));
  }

  void fill(T val) { std::fill(data_.begin(), data_.end(), val); }

  //--------------------------------------------------------------------------
  // View accessors

  //! Row i of a 2D+ tensor (contiguous 1D view)
  View1D<stored_type> row(size_t i)
  {
    auto cols = shape_[shape_.size() - 1];
    return {data_.data() + i * cols, cols, 1};
  }
  View1D<const stored_type> row(size_t i) const
  {
    auto cols = shape_[shape_.size() - 1];
    return {data_.data() + i * cols, cols, 1};
  }

  //! Column j of a 2D tensor (strided 1D view)
  View1D<stored_type> col(size_t j)
  {
    return {data_.data() + j, shape_[0], shape_[1]};
  }
  View1D<const stored_type> col(size_t j) const
  {
    return {data_.data() + j, shape_[0], shape_[1]};
  }

  //! Subrange of a 1D tensor
  View1D<stored_type> slice(size_t start, size_t end)
  {
    return {data_.data() + start, end - start, 1};
  }
  View1D<const stored_type> slice(size_t start, size_t end) const
  {
    return {data_.data() + start, end - start, 1};
  }

  //! Subrange to end of a 1D tensor
  View1D<stored_type> slice(size_t start)
  {
    return {data_.data() + start, data_.size() - start, 1};
  }
  View1D<const stored_type> slice(size_t start) const
  {
    return {data_.data() + start, data_.size() - start, 1};
  }

  //! Flat 1D view of all elements
  ViewFlat<stored_type> flat()
  {
    return {data_.data(), data_.size()};
  }
  ViewFlat<const stored_type> flat() const
  {
    return {data_.data(), data_.size()};
  }

  //--------------------------------------------------------------------------
  // Reductions

  T sum() const
  {
    T s = T(0);
    for (size_t i = 0; i < data_.size(); ++i)
      s += data_[i];
    return s;
  }

  //! Sum along an axis, reducing rank by 1 (defined out-of-line below)
  Tensor<T> sum(size_t axis) const;

  T prod() const
  {
    T p = T(1);
    for (size_t i = 0; i < data_.size(); ++i)
      p *= data_[i];
    return p;
  }

  bool any() const
  {
    for (size_t i = 0; i < data_.size(); ++i)
      if (data_[i])
        return true;
    return false;
  }

  bool all() const
  {
    for (size_t i = 0; i < data_.size(); ++i)
      if (!data_[i])
        return false;
    return true;
  }

  size_t argmin() const
  {
    return static_cast<size_t>(
      std::distance(data_.data(),
        std::min_element(data_.data(), data_.data() + data_.size())));
  }

  //--------------------------------------------------------------------------
  // Flip

  Tensor flip(size_t axis) const
  {
    size_t outer_size = 1;
    for (size_t d = 0; d < axis; ++d)
      outer_size *= shape_[d];
    size_t axis_size = shape_[axis];
    size_t inner_size = 1;
    for (size_t d = axis + 1; d < shape_.size(); ++d)
      inner_size *= shape_[d];

    Tensor r(shape_);
    for (size_t o = 0; o < outer_size; ++o)
      for (size_t a = 0; a < axis_size; ++a)
        for (size_t i = 0; i < inner_size; ++i)
          r.data_[(o * axis_size + (axis_size - 1 - a)) * inner_size + i] =
            data_[(o * axis_size + a) * inner_size + i];
    return r;
  }

  //--------------------------------------------------------------------------
  // Compound assignment operators (scalar)

  Tensor& operator+=(T val)
  {
    for (auto& x : data_)
      x += val;
    return *this;
  }
  Tensor& operator-=(T val)
  {
    for (auto& x : data_)
      x -= val;
    return *this;
  }
  Tensor& operator*=(T val)
  {
    for (auto& x : data_)
      x *= val;
    return *this;
  }
  Tensor& operator/=(T val)
  {
    for (auto& x : data_)
      x /= val;
    return *this;
  }

  //--------------------------------------------------------------------------
  // Compound assignment operators (tensor)

  Tensor& operator+=(const Tensor& o)
  {
    for (size_t i = 0; i < data_.size(); ++i)
      data_[i] += o.data_[i];
    return *this;
  }

  //--------------------------------------------------------------------------
  // Element-wise binary operators (tensor op tensor)

  Tensor operator+(const Tensor& o) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] + o.data_[i];
    return r;
  }
  Tensor operator-(const Tensor& o) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] - o.data_[i];
    return r;
  }
  Tensor operator/(const Tensor& o) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] / o.data_[i];
    return r;
  }

  //--------------------------------------------------------------------------
  // Element-wise binary operators (tensor op scalar)

  Tensor operator+(T val) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] + val;
    return r;
  }
  Tensor operator-(T val) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] - val;
    return r;
  }
  Tensor operator*(T val) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] * val;
    return r;
  }

  //--------------------------------------------------------------------------
  // Element-wise comparison operators (return Tensor<bool>)

  Tensor<bool> operator<=(T val) const
  {
    Tensor<bool> r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] <= val;
    return r;
  }
  Tensor<bool> operator<(T val) const
  {
    Tensor<bool> r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] < val;
    return r;
  }
  Tensor<bool> operator>=(T val) const
  {
    Tensor<bool> r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] >= val;
    return r;
  }
  Tensor<bool> operator>(T val) const
  {
    Tensor<bool> r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] > val;
    return r;
  }
  Tensor<bool> operator<(const Tensor& o) const
  {
    Tensor<bool> r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] < o.data_[i];
    return r;
  }

private:
  size_t compute_size() const
  {
    size_t s = 1;
    for (auto d : shape_)
      s *= d;
    return s;
  }

  //--------------------------------------------------------------------------
  // Data members

  vector<size_t> shape_;
  vector<storage_type<T>> data_;
};

//==============================================================================
// Free operators (scalar op tensor)
//==============================================================================

template<typename T,
  typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Tensor<T> operator*(T val, const Tensor<T>& arr)
{
  return arr * val;
}

template<typename T,
  typename = std::enable_if_t<std::is_arithmetic<T>::value>>
Tensor<T> operator+(T val, const Tensor<T>& arr)
{
  return arr + val;
}

// Mixed-type arithmetic: Tensor<T1> op Tensor<T2> -> Tensor<double>
template<typename T1, typename T2,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
Tensor<double> operator*(const Tensor<T1>& a, const Tensor<T2>& b)
{
  Tensor<double> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] =
      static_cast<double>(a.data()[i]) * static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
Tensor<double> operator/(const Tensor<T1>& a, const Tensor<T2>& b)
{
  Tensor<double> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] =
      static_cast<double>(a.data()[i]) / static_cast<double>(b.data()[i]);
  return r;
}

//==============================================================================
// View1D deferred method definitions (need Tensor to be complete)
//==============================================================================

template<typename T>
template<typename U>
View1D<T>& View1D<T>::operator=(const Tensor<U>& other)
{
  for (size_t i = 0; i < size_; ++i)
    data_[i * stride_] = static_cast<T>(other.data()[i]);
  return *this;
}

template<typename T>
template<typename U>
View1D<T>& View1D<T>::operator+=(const Tensor<U>& o)
{
  for (size_t i = 0; i < size_; ++i)
    data_[i * stride_] += o.data()[i];
  return *this;
}

//==============================================================================
// Tensor<T>::sum(axis) — reduces one dimension
//==============================================================================

template<typename T>
Tensor<T> Tensor<T>::sum(size_t axis) const
{
  // Build output shape (all dims except the summed axis)
  vector<size_t> out_shape;
  for (size_t d = 0; d < shape_.size(); ++d)
    if (d != axis)
      out_shape.push_back(shape_[d]);

  // Split dimensions into three zones: outer | axis | inner
  size_t outer_size = 1;
  for (size_t d = 0; d < axis; ++d)
    outer_size *= shape_[d];
  size_t axis_size = shape_[axis];
  size_t inner_size = 1;
  for (size_t d = axis + 1; d < shape_.size(); ++d)
    inner_size *= shape_[d];

  Tensor<T> result(out_shape, T(0));
  for (size_t o = 0; o < outer_size; ++o)
    for (size_t a = 0; a < axis_size; ++a)
      for (size_t i = 0; i < inner_size; ++i)
        result.data()[o * inner_size + i] +=
          data_[(o * axis_size + a) * inner_size + i];

  return result;
}

//==============================================================================
// Fixed2D<T, R, C>: compile-time fixed 2D tensor.
//==============================================================================

template<typename T, size_t R, size_t C>
class Fixed2D {
public:
  using value_type = T;

  template<typename I0, typename I1>
  T& operator()(I0 i, I1 j)
  {
    return data_[static_cast<size_t>(i) * C +
                 static_cast<size_t>(j)];
  }
  template<typename I0, typename I1>
  const T& operator()(I0 i, I1 j) const
  {
    return data_[static_cast<size_t>(i) * C +
                 static_cast<size_t>(j)];
  }

  T* data() { return data_; }
  const T* data() const { return data_; }
  constexpr size_t size() const { return R * C; }
  std::array<size_t, 2> shape() const { return {R, C}; }

  void fill(T val) { std::fill(data_, data_ + R * C, val); }

  T* begin() { return data_; }
  T* end() { return data_ + R * C; }
  const T* begin() const { return data_; }
  const T* end() const { return data_ + R * C; }

  //! Column view
  View1D<T> col(size_t j) { return {data_ + j, R, C}; }
  View1D<const T> col(size_t j) const { return {data_ + j, R, C}; }

  //! Flat view
  ViewFlat<T> flat() { return {data_, R * C}; }
  ViewFlat<const T> flat() const { return {data_, R * C}; }

private:
  T data_[R * C] = {};
};

//==============================================================================
// Free functions
//==============================================================================

// zeros
template<typename T>
Tensor<T> zeros(std::initializer_list<size_t> shape)
{
  vector<size_t> s(shape);
  return Tensor<T>(std::move(s), T(0));
}

template<typename T>
Tensor<T> zeros(const vector<size_t>& shape)
{
  return Tensor<T>(shape, T(0));
}

// zeros_like
template<typename T>
Tensor<T> zeros_like(const Tensor<T>& o)
{
  return Tensor<T>(o.shape(), T(0));
}

// full_like
template<typename T, typename V>
Tensor<T> full_like(const Tensor<T>& o, V val)
{
  return Tensor<T>(o.shape(), static_cast<T>(val));
}

// linspace
template<typename T>
Tensor<T> linspace(T start, T stop, size_t n)
{
  Tensor<T> result({n});
  if (n < 2) {
    result[0] = start;
    return result;
  }
  for (size_t i = 0; i < n; ++i) {
    result[i] = start + static_cast<T>(i) * (stop - start) /
                          static_cast<T>(n - 1);
  }
  return result;
}

// concatenate (two 1D tensors)
template<typename T>
Tensor<T> concatenate(const Tensor<T>& a, const Tensor<T>& b)
{
  size_t total = a.size() + b.size();
  Tensor<T> result({total});
  std::copy(a.data(), a.data() + a.size(), result.data());
  std::copy(b.data(), b.data() + b.size(), result.data() + a.size());
  return result;
}

// Element-wise math
template<typename T>
Tensor<T> log(const Tensor<T>& a)
{
  Tensor<T> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::log(a.data()[i]);
  return r;
}

template<typename T>
Tensor<T> abs(const Tensor<T>& a)
{
  Tensor<T> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::abs(a.data()[i]);
  return r;
}

// where with tensor true_val and scalar false_val
template<typename T, typename V>
Tensor<T> where(
  const Tensor<bool>& cond, const Tensor<T>& true_val, V false_val)
{
  Tensor<T> r(cond.shape());
  for (size_t i = 0; i < cond.size(); ++i)
    r.data()[i] = cond.data()[i] ? true_val.data()[i]
                                 : static_cast<T>(false_val);
  return r;
}

// nan_to_num
template<typename T>
Tensor<T> nan_to_num(const Tensor<T>& a, T nan_val = T(0),
  T posinf_val = std::numeric_limits<T>::max(),
  T neginf_val = std::numeric_limits<T>::lowest())
{
  Tensor<T> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i) {
    T val = a.data()[i];
    if (std::isnan(val))
      r.data()[i] = nan_val;
    else if (std::isinf(val))
      r.data()[i] = val > 0 ? posinf_val : neginf_val;
    else
      r.data()[i] = val;
  }
  return r;
}

//==============================================================================
// Type traits
//==============================================================================

//! Type trait that is true for Tensor and Fixed2D.
//! Used by hdf5_interface.h to select the correct write_dataset overload.
template<typename T>
struct is_tensor : std::false_type {};

template<typename T>
struct is_tensor<Tensor<T>> : std::true_type {};

template<typename T, size_t R, size_t C>
struct is_tensor<Fixed2D<T, R, C>> : std::true_type {};

} // namespace tensor
} // namespace openmc

#endif // OPENMC_TENSOR_H
