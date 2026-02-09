#ifndef OPENMC_TENSOR_H
#define OPENMC_TENSOR_H

// Drop-in replacement for xtensor functionality used by OpenMC.
// Only implements the subset of xtensor that OpenMC actually calls.
// Built on openmc::vector for future GPU portability.

#include "openmc/vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <tuple>
#include <type_traits>

namespace xt {

//==============================================================================
// Forward declarations
//==============================================================================

template<typename T, std::size_t N>
class xtensor;

template<typename T>
class xarray;

template<typename T, typename Shape>
class xtensor_fixed;

template<std::size_t... Dims>
struct xshape {};

//==============================================================================
// Slice/view helper types
//==============================================================================

struct xall_type {};

inline xall_type all()
{
  return {};
}

struct xrange_type {
  std::ptrdiff_t start;
  std::ptrdiff_t stop;
};

inline xrange_type range(std::ptrdiff_t start, std::ptrdiff_t stop)
{
  return {start, stop};
}

namespace placeholders {
constexpr std::ptrdiff_t _ = std::numeric_limits<std::ptrdiff_t>::max();
}

// Ownership tags for adapt()
struct no_ownership {};
struct acquire_ownership {};

//==============================================================================
// xshape traits for xtensor_fixed
//==============================================================================

namespace detail {
template<typename Shape>
struct xshape_traits;

template<std::size_t D0, std::size_t D1>
struct xshape_traits<xshape<D0, D1>> {
  static constexpr std::size_t ndim = 2;
  static constexpr std::size_t total = D0 * D1;
  static constexpr std::size_t dim0 = D0;
  static constexpr std::size_t dim1 = D1;
};
} // namespace detail

//==============================================================================
// Storage type: avoids std::vector<bool> specialization
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
// 1D view: a strided view into a contiguous buffer
//==============================================================================

template<typename T>
class xtensor_view_1d {
  T* data_;
  std::size_t size_;
  std::size_t stride_;

public:
  using value_type = std::remove_const_t<T>;

  xtensor_view_1d(T* data, std::size_t size, std::size_t stride = 1)
    : data_(data), size_(size), stride_(stride)
  {}

  T& operator()(std::size_t i) { return data_[i * stride_]; }
  const T& operator()(std::size_t i) const { return data_[i * stride_]; }
  T& operator[](std::size_t i) { return data_[i * stride_]; }
  const T& operator[](std::size_t i) const { return data_[i * stride_]; }

  std::size_t size() const { return size_; }

  // Assignment from another 1D view
  template<typename U>
  xtensor_view_1d& operator=(const xtensor_view_1d<U>& other)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] = other(i);
    return *this;
  }

  // Assignment from xarray
  template<typename U>
  xtensor_view_1d& operator=(const xarray<U>& other)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] = static_cast<T>(other.data()[i]);
    return *this;
  }

  // Assignment from xtensor (any dimension, uses linear indexing)
  template<typename U, std::size_t M>
  xtensor_view_1d& operator=(const xtensor<U, M>& other)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] = other.data()[i];
    return *this;
  }

  // Assignment from scalar
  template<typename U>
  auto operator=(U val) ->
    std::enable_if_t<std::is_arithmetic<U>::value, xtensor_view_1d&>
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] = val;
    return *this;
  }

  // Compound assignment operators
  template<typename U, std::size_t M>
  xtensor_view_1d& operator+=(const xtensor<U, M>& o)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] += o.data()[i];
    return *this;
  }

  template<typename U>
  xtensor_view_1d& operator+=(const xtensor_view_1d<U>& o)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] += o(i);
    return *this;
  }

  xtensor_view_1d& operator/=(value_type val)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] /= val;
    return *this;
  }

  xtensor_view_1d& operator*=(value_type val)
  {
    for (std::size_t i = 0; i < size_; ++i)
      data_[i * stride_] *= val;
    return *this;
  }

  // Strided iterator
  class const_iterator {
    const T* ptr_;
    std::size_t stride_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_const_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    const_iterator(const T* ptr, std::size_t stride)
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
    friend const_iterator operator+(difference_type n, const const_iterator& it)
    {
      return it + n;
    }
  };

  class iterator {
    T* ptr_;
    std::size_t stride_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_const_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    iterator(T* ptr, std::size_t stride) : ptr_(ptr), stride_(stride) {}
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
    bool operator==(const iterator& other) const { return ptr_ == other.ptr_; }
    bool operator!=(const iterator& other) const { return ptr_ != other.ptr_; }
    bool operator<(const iterator& other) const { return ptr_ < other.ptr_; }
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

  T* data() { return data_; }
  const T* data() const { return data_; }
  std::size_t stride() const { return stride_; }
};

//==============================================================================
// 2D view: a strided 2D view into a contiguous buffer
//==============================================================================

template<typename T>
class xtensor_view_2d {
  T* data_;
  std::size_t shape0_, shape1_;
  std::size_t stride0_, stride1_;

public:
  using value_type = std::remove_const_t<T>;

  xtensor_view_2d(T* data, std::size_t s0, std::size_t s1, std::size_t st0,
    std::size_t st1)
    : data_(data), shape0_(s0), shape1_(s1), stride0_(st0), stride1_(st1)
  {}

  T& operator()(std::size_t i, std::size_t j)
  {
    return data_[i * stride0_ + j * stride1_];
  }
  const T& operator()(std::size_t i, std::size_t j) const
  {
    return data_[i * stride0_ + j * stride1_];
  }

  std::size_t size() const { return shape0_ * shape1_; }
  std::array<std::size_t, 2> shape() const { return {shape0_, shape1_}; }

  // Assignment from 2D tensor
  template<typename U, std::size_t M>
  auto operator=(const xtensor<U, M>& other) ->
    std::enable_if_t<M == 2, xtensor_view_2d&>
  {
    for (std::size_t i = 0; i < shape0_; ++i)
      for (std::size_t j = 0; j < shape1_; ++j)
        data_[i * stride0_ + j * stride1_] = other(i, j);
    return *this;
  }

  // Assignment from scalar
  template<typename U>
  auto operator=(U val) ->
    std::enable_if_t<std::is_arithmetic<U>::value, xtensor_view_2d&>
  {
    for (std::size_t i = 0; i < shape0_; ++i)
      for (std::size_t j = 0; j < shape1_; ++j)
        data_[i * stride0_ + j * stride1_] = val;
    return *this;
  }
};

//==============================================================================
// Flat view: wraps all elements for bulk assignment
//==============================================================================

template<typename T>
class xtensor_view_flat {
  T* data_;
  std::size_t size_;

public:
  xtensor_view_flat(T* data, std::size_t size) : data_(data), size_(size) {}

  template<typename U>
  auto operator=(U val) ->
    std::enable_if_t<std::is_arithmetic<U>::value, xtensor_view_flat&>
  {
    std::fill(data_, data_ + size_, static_cast<T>(val));
    return *this;
  }

  template<typename Container>
  auto operator=(const Container& other) ->
    std::enable_if_t<!std::is_arithmetic<Container>::value, xtensor_view_flat&>
  {
    std::copy(other.data(), other.data() + size_, data_);
    return *this;
  }
};

//==============================================================================
// xtensor<T, N>: N-dimensional tensor backed by openmc::vector
//==============================================================================

template<typename T, std::size_t N>
class xtensor {
  openmc::vector<storage_type<T>> data_;
  std::array<std::size_t, N> shape_;

  std::size_t compute_size() const
  {
    std::size_t s = 1;
    for (std::size_t i = 0; i < N; ++i)
      s *= shape_[i];
    return s;
  }

public:
  using value_type = T;
  using stored_type = storage_type<T>;
  using iterator = typename openmc::vector<stored_type>::iterator;
  using const_iterator = typename openmc::vector<stored_type>::const_iterator;
  using shape_type = std::array<std::size_t, N>;

  // Default constructor
  xtensor() { shape_.fill(0); }

  // Shape constructor (from initializer list of size_t)
  xtensor(std::initializer_list<std::size_t> shape)
  {
    std::copy(shape.begin(), shape.end(), shape_.begin());
    data_.resize(compute_size());
  }

  // Shape + fill value constructor (from initializer list)
  xtensor(std::initializer_list<std::size_t> shape, T val)
  {
    std::copy(shape.begin(), shape.end(), shape_.begin());
    data_.assign(compute_size(), val);
  }

  // Shape constructor (from array)
  explicit xtensor(const std::array<std::size_t, N>& shape) : shape_(shape)
  {
    data_.resize(compute_size());
  }

  // Shape + fill value constructor (from array)
  xtensor(const std::array<std::size_t, N>& shape, T val) : shape_(shape)
  {
    data_.assign(compute_size(), val);
  }

  // Construct from data + shape
  xtensor(
    openmc::vector<stored_type>&& data, const std::array<std::size_t, N>& shape)
    : data_(std::move(data)), shape_(shape)
  {}

  xtensor(const openmc::vector<stored_type>& data,
    const std::array<std::size_t, N>& shape)
    : data_(data), shape_(shape)
  {}

  // 1D initializer_list constructor (only for N==1, T != size_t to avoid
  // ambiguity)
  template<std::size_t M = N, typename TT = T,
    typename = std::enable_if_t<M == 1 && !std::is_same<TT, std::size_t>::value>>
  xtensor(std::initializer_list<T> vals)
  {
    shape_[0] = vals.size();
    data_.assign(vals.begin(), vals.end());
  }

  // Static factory: from_shape
  static xtensor from_shape(const std::array<std::size_t, N>& shape)
  {
    return xtensor(shape);
  }

  static xtensor from_shape(std::initializer_list<std::size_t> shape)
  {
    xtensor result;
    std::copy(shape.begin(), shape.end(), result.shape_.begin());
    result.data_.resize(result.compute_size());
    return result;
  }

  // Multi-dimensional indexing
  template<typename... Indices>
  stored_type& operator()(Indices... indices)
  {
    return data_[offset(static_cast<std::size_t>(indices)...)];
  }

  template<typename... Indices>
  const stored_type& operator()(Indices... indices) const
  {
    return data_[offset(static_cast<std::size_t>(indices)...)];
  }

  // Linear indexing
  stored_type& operator[](std::size_t i) { return data_[i]; }
  const stored_type& operator[](std::size_t i) const { return data_[i]; }

  // Data access
  stored_type* data() { return data_.data(); }
  const stored_type* data() const { return data_.data(); }
  std::size_t size() const { return data_.size(); }
  const shape_type& shape() const { return shape_; }
  std::size_t shape(std::size_t dim) const { return shape_[dim]; }
  bool empty() const { return data_.empty(); }
  static constexpr std::size_t dimension() { return N; }

  // Resize
  void resize(std::initializer_list<std::size_t> shape)
  {
    std::copy(shape.begin(), shape.end(), shape_.begin());
    data_.resize(compute_size());
  }

  void resize(const std::array<std::size_t, N>& shape)
  {
    shape_ = shape;
    data_.resize(compute_size());
  }

  // Resize from vector<size_t> (needed for xarray compatibility)
  void resize(const openmc::vector<std::size_t>& shape)
  {
    for (std::size_t i = 0; i < N && i < shape.size(); ++i)
      shape_[i] = shape[i];
    data_.resize(compute_size());
  }

  // Resize from vector<hsize_t>
  void resize(const openmc::vector<unsigned long long>& shape)
  {
    for (std::size_t i = 0; i < N && i < shape.size(); ++i)
      shape_[i] = static_cast<std::size_t>(shape[i]);
    data_.resize(compute_size());
  }

  // Fill
  void fill(T val) { std::fill(data_.begin(), data_.end(), val); }

  // Iterators
  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }
  const_iterator cbegin() const { return data_.cbegin(); }
  const_iterator cend() const { return data_.cend(); }

  // Compound assignment with scalar
  xtensor& operator+=(T val)
  {
    for (auto& x : data_)
      x += val;
    return *this;
  }
  xtensor& operator-=(T val)
  {
    for (auto& x : data_)
      x -= val;
    return *this;
  }
  xtensor& operator*=(T val)
  {
    for (auto& x : data_)
      x *= val;
    return *this;
  }
  xtensor& operator/=(T val)
  {
    for (auto& x : data_)
      x /= val;
    return *this;
  }

  // Compound assignment with tensor
  xtensor& operator+=(const xtensor& o)
  {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] += o.data_[i];
    return *this;
  }
  xtensor& operator-=(const xtensor& o)
  {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] -= o.data_[i];
    return *this;
  }
  xtensor& operator*=(const xtensor& o)
  {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] *= o.data_[i];
    return *this;
  }
  xtensor& operator/=(const xtensor& o)
  {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] /= o.data_[i];
    return *this;
  }

  // Compound assignment with 1D view
  template<typename U>
  xtensor& operator+=(const xtensor_view_1d<U>& o)
  {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] += o(i);
    return *this;
  }

  // Binary ops: tensor op tensor (same type)
  xtensor operator+(const xtensor& o) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] + o.data_[i];
    return r;
  }
  xtensor operator-(const xtensor& o) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] - o.data_[i];
    return r;
  }
  xtensor operator*(const xtensor& o) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] * o.data_[i];
    return r;
  }
  xtensor operator/(const xtensor& o) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] / o.data_[i];
    return r;
  }

  // Binary ops: tensor op scalar
  xtensor operator+(T val) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] + val;
    return r;
  }
  xtensor operator-(T val) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] - val;
    return r;
  }
  xtensor operator*(T val) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] * val;
    return r;
  }
  xtensor operator/(T val) const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] / val;
    return r;
  }

  // Unary minus
  xtensor operator-() const
  {
    xtensor r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = -data_[i];
    return r;
  }

  // Comparison operators (element-wise, return bool tensor)
  xtensor<bool, N> operator<=(T val) const
  {
    xtensor<bool, N> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] <= val;
    return r;
  }
  xtensor<bool, N> operator<(T val) const
  {
    xtensor<bool, N> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] < val;
    return r;
  }
  xtensor<bool, N> operator>=(T val) const
  {
    xtensor<bool, N> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] >= val;
    return r;
  }
  xtensor<bool, N> operator>(T val) const
  {
    xtensor<bool, N> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] > val;
    return r;
  }

  // Tensor-tensor comparison
  xtensor<bool, N> operator<(const xtensor& o) const
  {
    xtensor<bool, N> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] < o.data_[i];
    return r;
  }

  // Conversion from 1D view
  template<typename U, std::size_t M = N,
    typename = std::enable_if_t<M == 1>>
  xtensor(const xtensor_view_1d<U>& v)
  {
    shape_[0] = v.size();
    data_.resize(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
      data_[i] = v(i);
  }

  // Conversion from 2D view
  template<typename U, std::size_t M = N,
    typename = std::enable_if_t<M == 2>>
  xtensor(const xtensor_view_2d<U>& v)
  {
    auto s = v.shape();
    shape_[0] = s[0];
    shape_[1] = s[1];
    data_.resize(s[0] * s[1]);
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        data_[i * s[1] + j] = v(i, j);
  }

  // Converting constructor from xarray
  xtensor(const xarray<T>& other);

  // Assignment from xarray
  xtensor& operator=(const xarray<T>& other);

  // Cross-type assignment from xtensor<U, N>
  template<typename U,
    typename = std::enable_if_t<!std::is_same<U, T>::value>>
  xtensor& operator=(const xtensor<U, N>& other)
  {
    shape_ = other.shape();
    data_.resize(other.size());
    for (std::size_t i = 0; i < other.size(); ++i)
      data_[i] = static_cast<T>(other.data()[i]);
    return *this;
  }

  // Assignment from 1D view
  template<typename U, std::size_t M = N, typename = std::enable_if_t<M == 1>>
  xtensor& operator=(const xtensor_view_1d<U>& v)
  {
    shape_[0] = v.size();
    data_.resize(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
      data_[i] = v(i);
    return *this;
  }

  // Assignment from initializer list of values (for 1D only)
  template<std::size_t M = N, typename TT = T,
    typename = std::enable_if_t<M == 1 && !std::is_same<TT, std::size_t>::value>>
  xtensor& operator=(std::initializer_list<T> vals)
  {
    shape_[0] = vals.size();
    data_.assign(vals.begin(), vals.end());
    return *this;
  }

private:
  // Offset calculations for row-major layout
  std::size_t offset(std::size_t i0) const { return i0; }

  std::size_t offset(std::size_t i0, std::size_t i1) const
  {
    return i0 * shape_[1] + i1;
  }

  std::size_t offset(std::size_t i0, std::size_t i1, std::size_t i2) const
  {
    return (i0 * shape_[1] + i1) * shape_[2] + i2;
  }

  std::size_t offset(
    std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return ((i0 * shape_[1] + i1) * shape_[2] + i2) * shape_[3] + i3;
  }
};

// scalar op tensor (same type)
template<typename T, std::size_t N,
  typename = std::enable_if_t<std::is_arithmetic<T>::value>>
xtensor<T, N> operator*(T val, const xtensor<T, N>& arr)
{
  return arr * val;
}

template<typename T, std::size_t N,
  typename = std::enable_if_t<std::is_arithmetic<T>::value>>
xtensor<T, N> operator+(T val, const xtensor<T, N>& arr)
{
  return arr + val;
}

template<typename T, std::size_t N>
xtensor<T, N> operator-(T val, const xtensor<T, N>& arr)
{
  xtensor<T, N> r(arr.shape());
  for (std::size_t i = 0; i < arr.size(); ++i)
    r.data()[i] = val - arr.data()[i];
  return r;
}

template<typename T, std::size_t N>
xtensor<T, N> operator/(T val, const xtensor<T, N>& arr)
{
  xtensor<T, N> r(arr.shape());
  for (std::size_t i = 0; i < arr.size(); ++i)
    r.data()[i] = val / arr.data()[i];
  return r;
}

// Mixed-type arithmetic: xtensor<T1> op xtensor<T2>
// Returns xtensor<double, N> for common cases (int*double, double/int, etc.)
template<typename T1, typename T2, std::size_t N,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
xtensor<double, N> operator*(
  const xtensor<T1, N>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) *
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
xtensor<double, N> operator/(
  const xtensor<T1, N>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) /
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
xtensor<double, N> operator+(
  const xtensor<T1, N>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) +
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N,
  typename = std::enable_if_t<!std::is_same<T1, T2>::value>>
xtensor<double, N> operator-(
  const xtensor<T1, N>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) -
                  static_cast<double>(b.data()[i]);
  return r;
}

//==============================================================================
// xarray<T>: dynamic-dimension tensor backed by openmc::vector
//==============================================================================

template<typename T>
class xarray {
  openmc::vector<storage_type<T>> data_;
  openmc::vector<std::size_t> shape_;

public:
  using value_type = T;
  using stored_type = storage_type<T>;
  using iterator = typename openmc::vector<stored_type>::iterator;
  using const_iterator = typename openmc::vector<stored_type>::const_iterator;

  xarray() = default;

  explicit xarray(const openmc::vector<std::size_t>& shape) : shape_(shape)
  {
    std::size_t total = 1;
    for (auto d : shape_)
      total *= d;
    data_.resize(total);
  }

  xarray(const openmc::vector<std::size_t>& shape, T val) : shape_(shape)
  {
    std::size_t total = 1;
    for (auto d : shape_)
      total *= d;
    data_.assign(total, val);
  }

  // Construct from initializer_list of values (creates 1D xarray)
  xarray(std::initializer_list<T> vals) : shape_({vals.size()})
  {
    data_.assign(vals.begin(), vals.end());
  }

  // Converting constructor from xtensor<T, N>
  template<std::size_t N>
  xarray(const xtensor<T, N>& t)
  {
    for (std::size_t i = 0; i < N; ++i)
      shape_.push_back(t.shape()[i]);
    data_.assign(t.data(), t.data() + t.size());
  }

  // Construct from vector<hsize_t> shape
  explicit xarray(const openmc::vector<unsigned long long>& shape)
  {
    for (auto d : shape)
      shape_.push_back(static_cast<std::size_t>(d));
    std::size_t total = 1;
    for (auto d : shape_)
      total *= d;
    data_.resize(total);
  }

  // Multi-dimensional indexing (up to 4D)
  template<typename... Indices>
  stored_type& operator()(Indices... indices)
  {
    return data_[compute_offset(static_cast<std::size_t>(indices)...)];
  }

  template<typename... Indices>
  const stored_type& operator()(Indices... indices) const
  {
    return data_[compute_offset(static_cast<std::size_t>(indices)...)];
  }

  stored_type& operator[](std::size_t i) { return data_[i]; }
  const stored_type& operator[](std::size_t i) const { return data_[i]; }

  stored_type* data() { return data_.data(); }
  const stored_type* data() const { return data_.data(); }
  std::size_t size() const { return data_.size(); }
  const openmc::vector<std::size_t>& shape() const { return shape_; }
  std::size_t shape(std::size_t dim) const { return shape_[dim]; }
  bool empty() const { return data_.empty(); }

  void resize(const openmc::vector<std::size_t>& shape)
  {
    shape_ = shape;
    std::size_t total = 1;
    for (auto d : shape_)
      total *= d;
    data_.resize(total);
  }

  void resize(const openmc::vector<unsigned long long>& shape)
  {
    shape_.clear();
    for (auto d : shape)
      shape_.push_back(static_cast<std::size_t>(d));
    std::size_t total = 1;
    for (auto d : shape_)
      total *= d;
    data_.resize(total);
  }

  // reshape: change shape without changing data (total must match)
  template<typename ShapeType>
  void reshape(const ShapeType& new_shape)
  {
    shape_.clear();
    for (auto d : new_shape)
      shape_.push_back(static_cast<std::size_t>(d));
  }

  void fill(T val) { std::fill(data_.begin(), data_.end(), val); }

  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }

  // Compound assignment
  xarray& operator+=(T val)
  {
    for (auto& x : data_)
      x += val;
    return *this;
  }
  xarray& operator-=(T val)
  {
    for (auto& x : data_)
      x -= val;
    return *this;
  }
  xarray& operator*=(T val)
  {
    for (auto& x : data_)
      x *= val;
    return *this;
  }
  xarray& operator/=(T val)
  {
    for (auto& x : data_)
      x /= val;
    return *this;
  }

  // Binary ops with scalar
  xarray operator*(T val) const
  {
    xarray r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] * val;
    return r;
  }
  xarray operator/(T val) const
  {
    xarray r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] / val;
    return r;
  }
  xarray operator+(T val) const
  {
    xarray r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] + val;
    return r;
  }
  xarray operator-(T val) const
  {
    xarray r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] - val;
    return r;
  }
  xarray operator-() const
  {
    xarray r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = -data_[i];
    return r;
  }

  // Comparison with scalar
  xarray<bool> operator>(T val) const
  {
    xarray<bool> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] > val;
    return r;
  }
  xarray<bool> operator<(T val) const
  {
    xarray<bool> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] < val;
    return r;
  }
  xarray<bool> operator<=(T val) const
  {
    xarray<bool> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] <= val;
    return r;
  }
  xarray<bool> operator>=(T val) const
  {
    xarray<bool> r(shape_);
    for (std::size_t i = 0; i < data_.size(); ++i)
      r.data()[i] = data_[i] >= val;
    return r;
  }

  // Convert to xtensor (explicit shape)
  template<std::size_t M>
  operator xtensor<T, M>() const
  {
    std::array<std::size_t, M> s {};
    for (std::size_t i = 0; i < M && i < shape_.size(); ++i)
      s[i] = shape_[i];
    xtensor<T, M> result(s);
    std::copy(data_.begin(), data_.end(), result.data());
    return result;
  }

private:
  std::size_t compute_offset(std::size_t i0) const { return i0; }

  std::size_t compute_offset(std::size_t i0, std::size_t i1) const
  {
    return i0 * shape_[1] + i1;
  }

  std::size_t compute_offset(
    std::size_t i0, std::size_t i1, std::size_t i2) const
  {
    return (i0 * shape_[1] + i1) * shape_[2] + i2;
  }

  std::size_t compute_offset(
    std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return ((i0 * shape_[1] + i1) * shape_[2] + i2) * shape_[3] + i3;
  }
};

// scalar op xarray
template<typename T>
xarray<T> operator*(T val, const xarray<T>& arr)
{
  return arr * val;
}

// xtensor converting constructor from xarray (defined after xarray is complete)
template<typename T, std::size_t N>
xtensor<T, N>::xtensor(const xarray<T>& other)
{
  for (std::size_t i = 0; i < N; ++i)
    shape_[i] = other.shape()[i];
  data_.assign(other.data(), other.data() + other.size());
}

// xtensor assignment from xarray (defined after xarray is complete)
template<typename T, std::size_t N>
xtensor<T, N>& xtensor<T, N>::operator=(const xarray<T>& other)
{
  for (std::size_t i = 0; i < N; ++i)
    shape_[i] = other.shape()[i];
  data_.assign(other.data(), other.data() + other.size());
  return *this;
}

// Mixed-type arithmetic: xarray<T1> op xtensor<T2, N>
// These handle cases like xarray<int> * xtensor<double, 1>
template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator*(const xarray<T1>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(b.shape());
  for (std::size_t i = 0; i < b.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) *
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator*(const xtensor<T1, N>& a, const xarray<T2>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) *
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator/(const xtensor<T1, N>& a, const xarray<T2>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) /
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator/(const xarray<T1>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(b.shape());
  for (std::size_t i = 0; i < b.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) /
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator+(const xarray<T1>& a, const xtensor<T2, N>& b)
{
  xtensor<double, N> r(b.shape());
  for (std::size_t i = 0; i < b.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) +
                  static_cast<double>(b.data()[i]);
  return r;
}

template<typename T1, typename T2, std::size_t N>
xtensor<double, N> operator+(const xtensor<T1, N>& a, const xarray<T2>& b)
{
  xtensor<double, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = static_cast<double>(a.data()[i]) +
                  static_cast<double>(b.data()[i]);
  return r;
}

//==============================================================================
// xtensor_fixed<T, xshape<Dims...>>: fixed-size tensor
//==============================================================================

template<typename T, typename Shape>
class xtensor_fixed {
  static constexpr std::size_t total_ = detail::xshape_traits<Shape>::total;
  static constexpr std::size_t dim0_ = detail::xshape_traits<Shape>::dim0;
  static constexpr std::size_t dim1_ = detail::xshape_traits<Shape>::dim1;

  T data_[total_] = {};

public:
  using value_type = T;

  template<typename I0, typename I1>
  T& operator()(I0 i, I1 j)
  {
    return data_[static_cast<std::size_t>(i) * dim1_ +
                 static_cast<std::size_t>(j)];
  }
  template<typename I0, typename I1>
  const T& operator()(I0 i, I1 j) const
  {
    return data_[static_cast<std::size_t>(i) * dim1_ +
                 static_cast<std::size_t>(j)];
  }

  T* data() { return data_; }
  const T* data() const { return data_; }
  std::size_t size() const { return total_; }
  std::array<std::size_t, 2> shape() const { return {dim0_, dim1_}; }
  void fill(T val) { std::fill(data_, data_ + total_, val); }
};

//==============================================================================
// adapt() functions
//==============================================================================

// Adapt a std::vector into a 1D xtensor (copy)
template<typename T>
xtensor<T, 1> adapt(const std::vector<T>& vec)
{
  xtensor<T, 1> result({vec.size()});
  std::copy(vec.begin(), vec.end(), result.data());
  return result;
}

// Adapt a vector with explicit shape into xarray (copy)
template<typename T, typename ShapeType>
xarray<T> adapt(const openmc::vector<T>& vec, const ShapeType& shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  xarray<T> result(s);
  std::copy(vec.begin(), vec.end(), result.data());
  return result;
}

// Adapt std::array with explicit shape (copy)
template<typename T, std::size_t M, typename ShapeType>
xarray<T> adapt(const std::array<T, M>& arr, const ShapeType& shape)
{
  openmc::vector<std::size_t> s;
  std::size_t total = 1;
  for (auto d : shape) {
    auto dim = static_cast<std::size_t>(d);
    s.push_back(dim);
    total *= dim;
  }
  xarray<T> result(s);
  std::copy(arr.data(), arr.data() + std::min(total, M), result.data());
  return result;
}

// Adapt vector with initializer_list shape (for brace-init-list deduction)
template<typename T>
xarray<T> adapt(const openmc::vector<T>& vec, std::initializer_list<int> shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  xarray<T> result(s);
  std::copy(vec.begin(), vec.end(), result.data());
  return result;
}

// Adapt std::array with initializer_list shape (for brace-init-list deduction)
template<typename T, std::size_t M>
xarray<T> adapt(const std::array<T, M>& arr, std::initializer_list<int> shape)
{
  openmc::vector<std::size_t> s;
  std::size_t total = 1;
  for (auto d : shape) {
    auto dim = static_cast<std::size_t>(d);
    s.push_back(dim);
    total *= dim;
  }
  xarray<T> result(s);
  std::copy(arr.data(), arr.data() + std::min(total, M), result.data());
  return result;
}

// Adapt raw pointer with no_ownership into xtensor (copy)
template<typename T, typename ShapeType>
xarray<T> adapt(
  const T* ptr, std::size_t n, no_ownership, const ShapeType& shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  xarray<T> result(s);
  std::copy(ptr, ptr + n, result.data());
  return result;
}

// Adapt raw pointer with acquire_ownership into xarray (copy + delete)
template<typename T, typename ShapeType>
xarray<T> adapt(
  T* ptr, std::size_t total, acquire_ownership, const ShapeType& shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  xarray<T> result(s);
  std::copy(ptr, ptr + total, result.data());
  delete[] ptr;
  return result;
}

// Adapt const pointer with shape into xarray (copy)
template<typename T, typename ShapeType>
xarray<T> adapt(const T* ptr, const ShapeType& shape)
{
  openmc::vector<std::size_t> s;
  std::size_t total = 1;
  for (auto d : shape) {
    auto dim = static_cast<std::size_t>(d);
    s.push_back(dim);
    total *= dim;
  }
  xarray<T> result(s);
  std::copy(ptr, ptr + total, result.data());
  return result;
}

//==============================================================================
// Construction helpers
//==============================================================================

template<typename T>
xarray<T> zeros(std::initializer_list<std::size_t> shape)
{
  openmc::vector<std::size_t> s(shape);
  xarray<T> result(s, T(0));
  return result;
}

// zeros with vector<size_t> shape
template<typename T>
xarray<T> zeros(const openmc::vector<std::size_t>& shape)
{
  return xarray<T>(shape, T(0));
}

// zeros with vector of any int type
template<typename T, typename I,
  typename = std::enable_if_t<std::is_integral<I>::value &&
    !std::is_same<I, std::size_t>::value>>
xarray<T> zeros(const openmc::vector<I>& shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  return xarray<T>(s, T(0));
}

template<typename T>
xtensor<T, 2> zeros(const std::array<std::size_t, 2>& shape)
{
  return xtensor<T, 2>(shape, T(0));
}

template<typename T>
xtensor<T, 3> zeros(const std::array<std::size_t, 3>& shape)
{
  return xtensor<T, 3>(shape, T(0));
}

template<typename T>
xarray<T> empty(std::initializer_list<std::size_t> shape)
{
  openmc::vector<std::size_t> s(shape);
  return xarray<T>(s);
}

// empty with int initializer_list (avoids narrowing errors)
template<typename T>
xarray<T> empty(std::initializer_list<int> shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  return xarray<T>(s);
}

template<typename T>
xtensor<T, 2> empty(const std::array<std::size_t, 2>& shape)
{
  return xtensor<T, 2>(shape);
}

template<typename T>
xtensor<T, 3> empty(const std::array<std::size_t, 3>& shape)
{
  return xtensor<T, 3>(shape);
}

// empty with std::array<int, N> shape
template<typename T, std::size_t N>
xarray<T> empty(const std::array<int, N>& shape)
{
  openmc::vector<std::size_t> s;
  for (auto d : shape)
    s.push_back(static_cast<std::size_t>(d));
  return xarray<T>(s);
}

template<typename T, std::size_t N>
xtensor<T, N> zeros_like(const xtensor<T, N>& o)
{
  return xtensor<T, N>(o.shape(), T(0));
}

template<typename T, std::size_t N>
xtensor<T, N> empty_like(const xtensor<T, N>& o)
{
  return xtensor<T, N>(o.shape());
}

template<typename T, std::size_t N>
xtensor<T, N> ones_like(const xtensor<T, N>& o)
{
  return xtensor<T, N>(o.shape(), T(1));
}

// full_like: create tensor with same shape, filled with value
template<typename T, std::size_t N, typename V>
xtensor<T, N> full_like(const xtensor<T, N>& o, V val)
{
  return xtensor<T, N>(o.shape(), static_cast<T>(val));
}

template<typename T>
xarray<T> empty_like(const xarray<T>& o)
{
  return xarray<T>(o.shape());
}

template<typename T>
xtensor<T, 1> linspace(T start, T stop, std::size_t n)
{
  xtensor<T, 1> result({n});
  if (n < 2) {
    result[0] = start;
    return result;
  }
  for (std::size_t i = 0; i < n; ++i) {
    result[i] = start + static_cast<T>(i) * (stop - start) /
                          static_cast<T>(n - 1);
  }
  return result;
}

//==============================================================================
// view() functions
//==============================================================================

// Resolve a range endpoint: if == placeholders::_, use dim as the end
inline std::size_t resolve_end(std::ptrdiff_t val, std::size_t dim)
{
  if (val == placeholders::_)
    return dim;
  if (val < 0)
    return dim + val;
  return static_cast<std::size_t>(val);
}

// view(2D tensor, range, int) -> 1D view of column slice
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 2>& a, xrange_type r, std::size_t col)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start * a.shape()[1] + col, stop - start, a.shape()[1]};
}
template<typename T>
xtensor_view_1d<const T> view(
  const xtensor<T, 2>& a, xrange_type r, std::size_t col)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start * a.shape()[1] + col, stop - start, a.shape()[1]};
}

// view(2D tensor, all, int) -> 1D column view
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 2>& a, xall_type, std::size_t col)
{
  return {a.data() + col, a.shape()[0], a.shape()[1]};
}
template<typename T>
xtensor_view_1d<const T> view(
  const xtensor<T, 2>& a, xall_type, std::size_t col)
{
  return {a.data() + col, a.shape()[0], a.shape()[1]};
}

// view(2D tensor, int) -> 1D row view
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 2>& a, std::size_t row)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xtensor<T, 2>& a, std::size_t row)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}

// view(1D tensor, range) -> 1D view
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 1>& a, xrange_type r)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start, stop - start, 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xtensor<T, 1>& a, xrange_type r)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start, stop - start, 1};
}

// view(3D tensor, int, int, all) -> 1D view along depth dimension
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 3>& a, std::size_t i, std::size_t j, xall_type)
{
  auto depth = a.shape()[2];
  return {a.data() + (i * a.shape()[1] + j) * depth, depth, 1};
}
template<typename T>
xtensor_view_1d<const T> view(
  const xtensor<T, 3>& a, std::size_t i, std::size_t j, xall_type)
{
  auto depth = a.shape()[2];
  return {a.data() + (i * a.shape()[1] + j) * depth, depth, 1};
}

// view(3D tensor, all, all, int) -> 2D view
template<typename T>
xtensor_view_2d<T> view(xtensor<T, 3>& a, xall_type, xall_type, std::size_t k)
{
  return {a.data() + k, a.shape()[0], a.shape()[1],
    a.shape()[1] * a.shape()[2], a.shape()[2]};
}
template<typename T>
xtensor_view_2d<const T> view(
  const xtensor<T, 3>& a, xall_type, xall_type, std::size_t k)
{
  return {a.data() + k, a.shape()[0], a.shape()[1],
    a.shape()[1] * a.shape()[2], a.shape()[2]};
}

// view(3D tensor, range, all, int) -> 2D view (subset of rows)
template<typename T>
xtensor_view_2d<T> view(
  xtensor<T, 3>& a, xrange_type r, xall_type, std::size_t k)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start * a.shape()[1] * a.shape()[2] + k, stop - start,
    a.shape()[1], a.shape()[1] * a.shape()[2], a.shape()[2]};
}
template<typename T>
xtensor_view_2d<const T> view(
  const xtensor<T, 3>& a, xrange_type r, xall_type, std::size_t k)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start * a.shape()[1] * a.shape()[2] + k, stop - start,
    a.shape()[1], a.shape()[1] * a.shape()[2], a.shape()[2]};
}

// view(fixed_2D, all, int) -> 1D column view
template<typename T, typename Shape>
xtensor_view_1d<T> view(xtensor_fixed<T, Shape>& a, xall_type, std::size_t col)
{
  auto s = a.shape();
  return {a.data() + col, s[0], s[1]};
}
template<typename T, typename Shape>
xtensor_view_1d<const T> view(
  const xtensor_fixed<T, Shape>& a, xall_type, std::size_t col)
{
  auto s = a.shape();
  return {a.data() + col, s[0], s[1]};
}

// view(fixed, all) -> flat view (for bulk assignment)
template<typename T, typename Shape>
xtensor_view_flat<T> view(xtensor_fixed<T, Shape>& a, xall_type)
{
  return {a.data(), a.size()};
}

// view(tensor, all) -> flat view (for bulk assignment)
template<typename T, std::size_t N>
xtensor_view_flat<T> view(xtensor<T, N>& a, xall_type)
{
  return {a.data(), a.size()};
}

// view(xarray, int) -> row view of 2D xarray
template<typename T>
xtensor_view_1d<T> view(xarray<T>& a, std::size_t row)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xarray<T>& a, std::size_t row)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}

// view(xarray, range) -> 1D view for 1D xarray
template<typename T>
xtensor_view_1d<T> view(xarray<T>& a, xrange_type r)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start, stop - start, 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xarray<T>& a, xrange_type r)
{
  auto start = resolve_end(r.start, a.shape()[0]);
  auto stop = resolve_end(r.stop, a.shape()[0]);
  return {a.data() + start, stop - start, 1};
}

// view(2D tensor, row, xall_type) -> 1D row view (entire row)
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 2>& a, std::size_t row, xall_type)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xtensor<T, 2>& a, std::size_t row, xall_type)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}

// view(2D tensor, row, range) -> 1D view of row subset
template<typename T>
xtensor_view_1d<T> view(xtensor<T, 2>& a, std::size_t row, xrange_type r)
{
  auto cols = a.shape()[1];
  auto start = resolve_end(r.start, cols);
  auto stop = resolve_end(r.stop, cols);
  return {a.data() + row * cols + start, stop - start, 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xtensor<T, 2>& a, std::size_t row, xrange_type r)
{
  auto cols = a.shape()[1];
  auto start = resolve_end(r.start, cols);
  auto stop = resolve_end(r.stop, cols);
  return {a.data() + row * cols + start, stop - start, 1};
}

// view(2D xarray, row, xall_type) -> 1D row view (entire row)
template<typename T>
xtensor_view_1d<T> view(xarray<T>& a, std::size_t row, xall_type)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xarray<T>& a, std::size_t row, xall_type)
{
  return {a.data() + row * a.shape()[1], a.shape()[1], 1};
}

// view(2D xarray, row, range) -> 1D view of row subset
template<typename T>
xtensor_view_1d<T> view(xarray<T>& a, std::size_t row, xrange_type r)
{
  auto cols = a.shape()[1];
  auto start = resolve_end(r.start, cols);
  auto stop = resolve_end(r.stop, cols);
  return {a.data() + row * cols + start, stop - start, 1};
}
template<typename T>
xtensor_view_1d<const T> view(const xarray<T>& a, std::size_t row, xrange_type r)
{
  auto cols = a.shape()[1];
  auto start = resolve_end(r.start, cols);
  auto stop = resolve_end(r.stop, cols);
  return {a.data() + row * cols + start, stop - start, 1};
}

// row() and col() for 2D tensors
template<typename T>
xtensor_view_1d<T> row(xtensor<T, 2>& a, std::size_t i)
{
  return {a.data() + i * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<const T> row(const xtensor<T, 2>& a, std::size_t i)
{
  return {a.data() + i * a.shape()[1], a.shape()[1], 1};
}
template<typename T>
xtensor_view_1d<T> col(xtensor<T, 2>& a, std::size_t j)
{
  return {a.data() + j, a.shape()[0], a.shape()[1]};
}
template<typename T>
xtensor_view_1d<const T> col(const xtensor<T, 2>& a, std::size_t j)
{
  return {a.data() + j, a.shape()[0], a.shape()[1]};
}

//==============================================================================
// Math / reduction functions
//==============================================================================

// sum - returns a callable proxy that yields the scalar when called with ()
template<typename T>
struct sum_proxy {
  T value;
  T operator()() const { return value; }
  operator T() const { return value; }
};

template<typename T, std::size_t N>
sum_proxy<T> sum(const xtensor<T, N>& a)
{
  T s = T(0);
  for (std::size_t i = 0; i < a.size(); ++i)
    s += a.data()[i];
  return {s};
}

template<typename T>
sum_proxy<T> sum(const xtensor_view_1d<T>& v)
{
  T s = T(0);
  for (std::size_t i = 0; i < v.size(); ++i)
    s += v(i);
  return {s};
}

template<typename T>
sum_proxy<T> sum(const xarray<T>& a)
{
  T s = T(0);
  for (std::size_t i = 0; i < a.size(); ++i)
    s += a.data()[i];
  return {s};
}

// sum along axis - reduces one dimension
// 2D sum along axis -> 1D
template<typename T>
xtensor<T, 1> sum(const xtensor<T, 2>& a, std::initializer_list<int> axes)
{
  int axis = *axes.begin();
  auto s = a.shape();
  if (axis == 0) {
    xtensor<T, 1> r({s[1]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        r(j) += a(i, j);
    return r;
  } else {
    xtensor<T, 1> r({s[0]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        r(i) += a(i, j);
    return r;
  }
}

// 3D sum along axis -> 2D
template<typename T>
xtensor<T, 2> sum(const xtensor<T, 3>& a, std::initializer_list<int> axes)
{
  int axis = *axes.begin();
  auto s = a.shape();
  if (axis == 0) {
    xtensor<T, 2> r({s[1], s[2]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          r(j, k) += a(i, j, k);
    return r;
  } else if (axis == 1) {
    xtensor<T, 2> r({s[0], s[2]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          r(i, k) += a(i, j, k);
    return r;
  } else {
    xtensor<T, 2> r({s[0], s[1]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          r(i, j) += a(i, j, k);
    return r;
  }
}

// 4D sum along axis -> 3D
template<typename T>
xtensor<T, 3> sum(const xtensor<T, 4>& a, std::initializer_list<int> axes)
{
  int axis = *axes.begin();
  auto s = a.shape();
  if (axis == 3) {
    xtensor<T, 3> r({s[0], s[1], s[2]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          for (std::size_t l = 0; l < s[3]; ++l)
            r(i, j, k) += a(i, j, k, l);
    return r;
  } else if (axis == 2) {
    xtensor<T, 3> r({s[0], s[1], s[3]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          for (std::size_t l = 0; l < s[3]; ++l)
            r(i, j, l) += a(i, j, k, l);
    return r;
  } else if (axis == 1) {
    xtensor<T, 3> r({s[0], s[2], s[3]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          for (std::size_t l = 0; l < s[3]; ++l)
            r(i, k, l) += a(i, j, k, l);
    return r;
  } else {
    xtensor<T, 3> r({s[1], s[2], s[3]}, T(0));
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        for (std::size_t k = 0; k < s[2]; ++k)
          for (std::size_t l = 0; l < s[3]; ++l)
            r(j, k, l) += a(i, j, k, l);
    return r;
  }
}

// prod - returns a callable proxy (like sum_proxy) for prod(...)() pattern
template<typename T>
struct prod_proxy {
  T value;
  T operator()() const { return value; }
  operator T() const { return value; }
};

template<typename T, std::size_t N>
prod_proxy<T> prod(const xtensor<T, N>& a)
{
  T p = T(1);
  for (std::size_t i = 0; i < a.size(); ++i)
    p *= a.data()[i];
  return {p};
}

template<typename T>
prod_proxy<T> prod(const xarray<T>& a)
{
  T p = T(1);
  for (std::size_t i = 0; i < a.size(); ++i)
    p *= a.data()[i];
  return {p};
}

// any (for bool tensors)
template<std::size_t N>
bool any(const xtensor<bool, N>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a.data()[i])
      return true;
  return false;
}

template<typename T, std::size_t N>
bool any(const xtensor<T, N>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a.data()[i])
      return true;
  return false;
}

// any/all for xarray
template<typename T>
bool any(const xarray<T>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a.data()[i])
      return true;
  return false;
}

template<typename T>
bool all(const xarray<T>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (!a.data()[i])
      return false;
  return true;
}

// all (for bool tensors)
template<std::size_t N>
bool all(const xtensor<bool, N>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (!a.data()[i])
      return false;
  return true;
}

template<typename T, std::size_t N>
bool all(const xtensor<T, N>& a)
{
  for (std::size_t i = 0; i < a.size(); ++i)
    if (!a.data()[i])
      return false;
  return true;
}

// argmin - returns a subscriptable proxy for argmin(...)[0] pattern
struct argmin_result {
  std::size_t value;
  std::size_t operator[](std::size_t) const { return value; }
  operator std::size_t() const { return value; }
  operator int() const { return static_cast<int>(value); }
};

template<typename T, std::size_t N>
argmin_result argmin(const xtensor<T, N>& a)
{
  return {static_cast<std::size_t>(std::distance(a.data(),
    std::min_element(a.data(), a.data() + a.size())))};
}

template<typename T>
argmin_result argmin(const xarray<T>& a)
{
  return {static_cast<std::size_t>(std::distance(a.data(),
    std::min_element(a.data(), a.data() + a.size())))};
}

// Element-wise math
template<typename T, std::size_t N>
xtensor<T, N> log(const xtensor<T, N>& a)
{
  xtensor<T, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::log(a.data()[i]);
  return r;
}

template<typename T>
xarray<T> log(const xarray<T>& a)
{
  xarray<T> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::log(a.data()[i]);
  return r;
}

template<typename T, std::size_t N>
xtensor<T, N> exp(const xtensor<T, N>& a)
{
  xtensor<T, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::exp(a.data()[i]);
  return r;
}

template<typename T, std::size_t N>
xtensor<T, N> abs(const xtensor<T, N>& a)
{
  xtensor<T, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::abs(a.data()[i]);
  return r;
}

template<typename T>
xarray<T> abs(const xarray<T>& a)
{
  xarray<T> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::abs(a.data()[i]);
  return r;
}

// where(condition, true_val, false_val) - tensor condition
template<typename T, std::size_t N, typename U, typename V>
xtensor<T, N> where(const xtensor<bool, N>& cond, U true_val, V false_val)
{
  xtensor<T, N> r(cond.shape());
  for (std::size_t i = 0; i < cond.size(); ++i)
    r.data()[i] = cond.data()[i] ? static_cast<T>(true_val)
                                 : static_cast<T>(false_val);
  return r;
}

// where with tensor true_val and scalar false_val
template<typename T, std::size_t N, typename V>
xtensor<T, N> where(
  const xtensor<bool, N>& cond, const xtensor<T, N>& true_val, V false_val)
{
  xtensor<T, N> r(cond.shape());
  for (std::size_t i = 0; i < cond.size(); ++i)
    r.data()[i] = cond.data()[i] ? true_val.data()[i]
                                 : static_cast<T>(false_val);
  return r;
}

// where with scalar true_val and tensor false_val
template<typename T, std::size_t N, typename U>
xtensor<T, N> where(
  const xtensor<bool, N>& cond, U true_val, const xtensor<T, N>& false_val)
{
  xtensor<T, N> r(cond.shape());
  for (std::size_t i = 0; i < cond.size(); ++i)
    r.data()[i] = cond.data()[i] ? static_cast<T>(true_val)
                                 : false_val.data()[i];
  return r;
}

// where with tensor true_val and tensor false_val
template<typename T, std::size_t N>
xtensor<T, N> where(const xtensor<bool, N>& cond,
  const xtensor<T, N>& true_val, const xtensor<T, N>& false_val)
{
  xtensor<T, N> r(cond.shape());
  for (std::size_t i = 0; i < cond.size(); ++i)
    r.data()[i] = cond.data()[i] ? true_val.data()[i] : false_val.data()[i];
  return r;
}

// nan_to_num
template<typename T, std::size_t N>
xtensor<T, N> nan_to_num(const xtensor<T, N>& a, T nan_val = T(0),
  T posinf_val = std::numeric_limits<T>::max(),
  T neginf_val = std::numeric_limits<T>::lowest())
{
  xtensor<T, N> r(a.shape());
  for (std::size_t i = 0; i < a.size(); ++i) {
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

// eval: returns a copy (materializes lazy expressions)
template<typename T, std::size_t N>
xtensor<T, N> eval(const xtensor<T, N>& a)
{
  return a;
}

template<typename T>
xarray<T> eval(const xarray<T>& a)
{
  return a;
}

//==============================================================================
// concatenate & xtuple
//==============================================================================

// Generic xtuple_holder that stores (pointer, size) pairs for any container type
template<typename T>
struct xtuple_holder {
  struct entry {
    const T* ptr;
    std::size_t size;
  };
  std::array<entry, 10> entries;
  std::size_t count;
};

// xtuple for 2 args - accepts any container with data() and size()
template<typename A, typename B>
auto xtuple(const A& a, const B& b)
  -> xtuple_holder<std::remove_const_t<std::remove_pointer_t<decltype(a.data())>>>
{
  using T = std::remove_const_t<std::remove_pointer_t<decltype(a.data())>>;
  xtuple_holder<T> h {};
  h.entries[0] = {a.data(), a.size()};
  h.entries[1] = {b.data(), b.size()};
  h.count = 2;
  return h;
}

// xtuple for 3 args
template<typename A, typename B, typename C>
auto xtuple(const A& a, const B& b, const C& c)
  -> xtuple_holder<std::remove_const_t<std::remove_pointer_t<decltype(a.data())>>>
{
  using T = std::remove_const_t<std::remove_pointer_t<decltype(a.data())>>;
  xtuple_holder<T> h {};
  h.entries[0] = {a.data(), a.size()};
  h.entries[1] = {b.data(), b.size()};
  h.entries[2] = {c.data(), c.size()};
  h.count = 3;
  return h;
}

template<typename T>
xtensor<T, 1> concatenate(const xtuple_holder<T>& tup)
{
  std::size_t total = 0;
  for (std::size_t i = 0; i < tup.count; ++i)
    total += tup.entries[i].size;

  xtensor<T, 1> result({total});
  std::size_t pos = 0;
  for (std::size_t i = 0; i < tup.count; ++i) {
    std::copy(tup.entries[i].ptr, tup.entries[i].ptr + tup.entries[i].size,
      result.data() + pos);
    pos += tup.entries[i].size;
  }
  return result;
}

// flip (reverse 1D tensor)
template<typename T>
xtensor<T, 1> flip(const xtensor<T, 1>& a)
{
  xtensor<T, 1> r({a.size()});
  for (std::size_t i = 0; i < a.size(); ++i)
    r.data()[i] = a.data()[a.size() - 1 - i];
  return r;
}

// flip 2D along axis 0
template<typename T>
xtensor<T, 2> flip(const xtensor<T, 2>& a, std::size_t axis = 0)
{
  auto s = a.shape();
  xtensor<T, 2> r(s);
  if (axis == 0) {
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        r(i, j) = a(s[0] - 1 - i, j);
  } else {
    for (std::size_t i = 0; i < s[0]; ++i)
      for (std::size_t j = 0; j < s[1]; ++j)
        r(i, j) = a(i, s[1] - 1 - j);
  }
  return r;
}

//==============================================================================
// Compatibility stubs
//==============================================================================

// is_xt_container trait for SFINAE in write_dataset
template<typename T>
struct is_xt_container : std::false_type {};

template<typename T, std::size_t N>
struct is_xt_container<xtensor<T, N>> : std::true_type {};

template<typename T>
struct is_xt_container<xarray<T>> : std::true_type {};

template<typename T, typename S>
struct is_xt_container<xtensor_fixed<T, S>> : std::true_type {};

// Compatibility types (unused but referenced in some template code)
template<typename CP, typename... Args>
using xbuffer_adaptor = openmc::vector<std::remove_pointer_t<std::remove_reference_t<CP>>>;

template<typename EC, std::size_t N, typename... Args>
using xtensor_adaptor = xtensor<typename EC::value_type, N>;

// noalias: no-op passthrough (xtensor optimization hint)
template<typename T>
T& noalias(T& x)
{
  return x;
}

} // namespace xt

#endif // OPENMC_TENSOR_H
