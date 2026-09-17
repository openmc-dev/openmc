//! \file tensor.h
//! \brief Multi-dimensional tensor types for OpenMC.
//!
//! Tensor<T> is the primary type: a dynamic-rank owning container that stores
//! elements contiguously in row-major order.  View<T> is a lightweight
//! non-owning reference into a Tensor's storage, returned by the slice()
//! method and flat().  StaticTensor2D<T,R,C> is a small stack-allocated 2D
//! array used only for simulation::global_tallies.
//!
//! Slicing follows numpy conventions: each axis takes an index (rank-reducing),
//! All (keep entire axis), or Range (keep sub-range).  For example,
//! arr.slice(0, all, range(2, 5)) is equivalent to numpy's arr[0, :, 2:5].
//!
//! View is declared before Tensor because Tensor's methods return View objects.

#ifndef OPENMC_TENSOR_H
#define OPENMC_TENSOR_H

#include "openmc/vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
class StaticTensor2D;

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
// Slice argument types
//
// Used with the variadic slice() method on Tensor, View, and StaticTensor2D.
// Each argument corresponds to one axis: a plain integer fixes that axis at
// a single index (rank-reducing), All keeps the entire axis, and Range keeps
// a sub-range.
//==============================================================================

//! Keep an entire axis (equivalent to numpy's ':' or xtensor's xt::all())
struct All {};
constexpr All all {};

//! Sub-range along an axis [start, end)
struct Range {
  size_t start;
  size_t end; // SIZE_MAX means "to end of axis"
};

//! Create a Range [start, end)
inline Range range(size_t start, size_t end)
{
  return {start, end};
}

//! Create a Range [0, end)
inline Range range(size_t end)
{
  return {0, end};
}

namespace detail {

//! Internal: normalized representation of a per-axis slice argument
struct SliceArg {
  enum Kind { INDEX, ALL, RANGE } kind;
  size_t start;
  size_t end;
};

inline SliceArg to_slice_arg(All)
{
  return {SliceArg::ALL, 0, 0};
}
inline SliceArg to_slice_arg(Range r)
{
  return {SliceArg::RANGE, r.start, r.end};
}

template<typename I>
inline
  typename std::enable_if<std::is_integral<I>::value || std::is_enum<I>::value,
    SliceArg>::type
  to_slice_arg(I i)
{
  return {SliceArg::INDEX, static_cast<size_t>(i), 0};
}

//! Result of a slice computation: pointer offset + new shape/strides
struct SliceResult {
  size_t ptr_offset;
  vector<size_t> shape;
  vector<size_t> strides;
};

//! Compute the result of applying slice arguments to shape/strides
template<typename First, typename... Rest>
SliceResult compute_slice(const vector<size_t>& shape,
  const vector<size_t>& strides, First first, Rest... rest)
{
  const size_t n = 1 + sizeof...(Rest);
  SliceArg args[1 + sizeof...(Rest)] = {
    to_slice_arg(first), to_slice_arg(rest)...};

  size_t offset = 0;
  vector<size_t> new_shape;
  vector<size_t> new_strides;

  for (size_t a = 0; a < n; ++a) {
    switch (args[a].kind) {
    case SliceArg::INDEX:
      offset += args[a].start * strides[a];
      break;
    case SliceArg::ALL:
      new_shape.push_back(shape[a]);
      new_strides.push_back(strides[a]);
      break;
    case SliceArg::RANGE: {
      offset += args[a].start * strides[a];
      size_t end = (args[a].end == SIZE_MAX) ? shape[a] : args[a].end;
      new_shape.push_back(end - args[a].start);
      new_strides.push_back(strides[a]);
      break;
    }
    }
  }

  // Trailing axes not covered by arguments are implicitly All.
  // This matches numpy: a[i] on a 2D array returns a 1D row.
  for (size_t a = n; a < shape.size(); ++a) {
    new_shape.push_back(shape[a]);
    new_strides.push_back(strides[a]);
  }

  return {offset, std::move(new_shape), std::move(new_strides)};
}

} // namespace detail

//==============================================================================
// View<T>: a non-owning N-dimensional view into a tensor's storage.
//
// Holds a base pointer, shape, and strides (in elements).  Supports arbitrary
// rank and multi-axis slicing via the variadic slice() method.
//==============================================================================

template<typename T>
class View {
public:
  //--------------------------------------------------------------------------
  // Constructors

  View(T* data, vector<size_t> shape, vector<size_t> strides)
    : data_(data), shape_(std::move(shape)), strides_(std::move(strides))
  {}

  // Explicitly default copy/move constructors (declaring copy assignment
  // below would otherwise suppress the implicit move constructor).
  View(const View&) = default;
  View(View&&) = default;

  //--------------------------------------------------------------------------
  // Indexing

  //! Multi-index element access (1D, 2D, 3D, ...)
  template<typename... Indices>
  T& operator()(Indices... indices)
  {
    const size_t idx[] = {static_cast<size_t>(indices)...};
    size_t off = 0;
    for (size_t d = 0; d < sizeof...(Indices); ++d)
      off += idx[d] * strides_[d];
    return data_[off];
  }

  template<typename... Indices>
  const T& operator()(Indices... indices) const
  {
    const size_t idx[] = {static_cast<size_t>(indices)...};
    size_t off = 0;
    for (size_t d = 0; d < sizeof...(Indices); ++d)
      off += idx[d] * strides_[d];
    return data_[off];
  }

  //! Flat logical index (row-major order)
  T& operator[](size_t i) { return data_[flat_to_offset(i)]; }
  const T& operator[](size_t i) const { return data_[flat_to_offset(i)]; }

  //--------------------------------------------------------------------------
  // Accessors

  size_t size() const
  {
    size_t s = 1;
    for (auto d : shape_)
      s *= d;
    return s;
  }
  size_t ndim() const { return shape_.size(); }
  size_t shape(size_t axis) const { return shape_[axis]; }
  const vector<size_t>& shape_vec() const { return shape_; }
  T* data() { return data_; }
  const T* data() const { return data_; }

  //--------------------------------------------------------------------------
  // View accessors

  //! Multi-axis slice.  Each argument corresponds to one axis and is either:
  //!   - an integer (fixes that axis, rank-reducing)
  //!   - All (keeps entire axis)
  //!   - Range (keeps sub-range along that axis)
  //! Example: v.slice(0, all, range(2, 5)) == numpy v[0, :, 2:5]
  template<typename First, typename... Rest>
  View<T> slice(First first, Rest... rest)
  {
    auto r = detail::compute_slice(shape_, strides_, first, rest...);
    return {data_ + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }

  template<typename First, typename... Rest>
  View<const T> slice(First first, Rest... rest) const
  {
    auto r = detail::compute_slice(shape_, strides_, first, rest...);
    return {data_ + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }

  //--------------------------------------------------------------------------
  // Assignment operators

  //! Copy assignment: element-wise deep copy (writes through data pointer).
  //! Without this, the compiler's implicit copy assignment just copies the
  //! View metadata (pointer, shape, strides) instead of the viewed data.
  View& operator=(const View& other)
  {
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      data_[flat_to_offset(i)] = other[i];
    return *this;
  }

  //! Fill all elements with a scalar
  View& operator=(T val)
  {
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      data_[flat_to_offset(i)] = val;
    return *this;
  }

  //! Assignment from initializer_list (for 1D views)
  View& operator=(std::initializer_list<T> vals)
  {
    auto it = vals.begin();
    for (size_t i = 0; i < size() && it != vals.end(); ++i, ++it)
      data_[flat_to_offset(i)] = *it;
    return *this;
  }

  //! Assignment from another View
  template<typename U>
  View& operator=(const View<U>& other)
  {
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      data_[flat_to_offset(i)] = other[i];
    return *this;
  }

  //! Assignment from Tensor (forward-declared, defined after Tensor)
  template<typename U>
  View& operator=(const Tensor<U>& other);

  //! Compound addition from Tensor (forward-declared, defined after Tensor)
  template<typename U>
  View& operator+=(const Tensor<U>& o);

  //! Compound multiply by scalar
  View& operator*=(T val)
  {
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      data_[flat_to_offset(i)] *= val;
    return *this;
  }

  //! Compound divide by scalar
  View& operator/=(T val)
  {
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      data_[flat_to_offset(i)] /= val;
    return *this;
  }

  //--------------------------------------------------------------------------
  // Reductions

  //! Sum of all elements
  T sum() const
  {
    // remove_const needed so accumulator is mutable when T is const-qualified
    std::remove_const_t<T> s = 0;
    size_t n = size();
    for (size_t i = 0; i < n; ++i)
      s += data_[flat_to_offset(i)];
    return s;
  }

  //--------------------------------------------------------------------------
  // Iterators
  //
  // Lightweight row-major iterator parameterized on pointer type (Ptr).
  // Stores a flat logical position and converts to a physical offset on each
  // dereference via divmod over shape/strides.  For contiguous 1D views (the
  // common case) the divmod chain reduces to a single multiply-by-1, which
  // the compiler optimizes away.
  //
  // view_iterator<T*>       = mutable iterator  (from non-const View)
  // view_iterator<const T*> = read-only iterator (from const View)

  template<typename Ptr>
  class view_iterator {
    Ptr base_;
    size_t count_;
    const size_t* shape_;
    const size_t* strides_;
    size_t ndim_;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = std::remove_const_t<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = Ptr;
    using reference = decltype(*std::declval<Ptr>());

    view_iterator(Ptr base, size_t count, const View* v)
      : base_(base), count_(count), shape_(v->shape_.data()),
        strides_(v->strides_.data()), ndim_(v->shape_.size())
    {}

    reference operator*() const { return base_[offset()]; }
    reference operator[](difference_type n) const
    {
      return base_[offset_of(count_ + n)];
    }
    view_iterator& operator++()
    {
      ++count_;
      return *this;
    }
    view_iterator operator++(int)
    {
      auto tmp = *this;
      ++count_;
      return tmp;
    }
    view_iterator& operator--()
    {
      --count_;
      return *this;
    }
    view_iterator operator+(difference_type n) const
    {
      auto tmp = *this;
      tmp.count_ += n;
      return tmp;
    }
    view_iterator operator-(difference_type n) const
    {
      auto tmp = *this;
      tmp.count_ -= n;
      return tmp;
    }
    difference_type operator-(const view_iterator& o) const
    {
      return static_cast<difference_type>(count_) -
             static_cast<difference_type>(o.count_);
    }
    view_iterator& operator+=(difference_type n)
    {
      count_ += n;
      return *this;
    }
    view_iterator& operator-=(difference_type n)
    {
      count_ -= n;
      return *this;
    }
    bool operator==(const view_iterator& o) const { return count_ == o.count_; }
    bool operator!=(const view_iterator& o) const { return count_ != o.count_; }
    bool operator<(const view_iterator& o) const { return count_ < o.count_; }
    bool operator>(const view_iterator& o) const { return count_ > o.count_; }
    bool operator<=(const view_iterator& o) const { return count_ <= o.count_; }
    bool operator>=(const view_iterator& o) const { return count_ >= o.count_; }
    friend view_iterator operator+(difference_type n, const view_iterator& it)
    {
      return it + n;
    }

  private:
    size_t offset() const { return offset_of(count_); }
    size_t offset_of(size_t flat) const
    {
      size_t off = 0;
      for (int d = static_cast<int>(ndim_) - 1; d >= 0; --d) {
        off += (flat % shape_[d]) * strides_[d];
        flat /= shape_[d];
      }
      return off;
    }
  };

  using iterator = view_iterator<T*>;
  using const_iterator = view_iterator<const T*>;

  iterator begin() { return {data_, 0, this}; }
  iterator end() { return {data_, size(), this}; }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }
  const_iterator cbegin() const { return {data_, 0, this}; }
  const_iterator cend() const { return {data_, size(), this}; }

private:
  //! Convert a logical flat index (row-major) to a physical element offset
  size_t flat_to_offset(size_t flat) const
  {
    size_t off = 0;
    for (int d = static_cast<int>(shape_.size()) - 1; d >= 0; --d) {
      off += (flat % shape_[d]) * strides_[d];
      flat /= shape_[d];
    }
    return off;
  }

  T* data_;
  vector<size_t> shape_;
  vector<size_t> strides_;
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

  //! Construct with shape (uninitialized for arithmetic types via vector
  //! resize)
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

  //! 1D copy from raw pointer + count
  Tensor(const T* ptr, size_t count) : shape_({count}), data_(ptr, ptr + count)
  {}

  //! Copy from View (preserves view's shape)
  template<typename U>
  explicit Tensor(const View<U>& v) : shape_(v.shape_vec())
  {
    size_t n = v.size();
    data_.resize(n);
    for (size_t i = 0; i < n; ++i)
      data_[i] = v[i];
  }

  //--------------------------------------------------------------------------
  // Assignment

  //! Assignment from View
  template<typename U>
  Tensor& operator=(const View<U>& v)
  {
    shape_ = v.shape_vec();
    size_t n = v.size();
    data_.resize(n);
    for (size_t i = 0; i < n; ++i)
      data_[i] = v[i];
    return *this;
  }

  //! Assignment from initializer_list of values (1D)
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
  size_t shape(size_t dim) const
  {
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

  //! First and last element
  stored_type& front() { return data_.front(); }
  const stored_type& front() const { return data_.front(); }
  stored_type& back() { return data_.back(); }
  const stored_type& back() const { return data_.back(); }

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

  void reshape(const vector<size_t>& new_shape) { shape_ = new_shape; }

  void fill(T val) { std::fill(data_.begin(), data_.end(), val); }

  //--------------------------------------------------------------------------
  // View accessors

  //! Fix one axis at a given index, returning an (N-1)-dimensional view
  //! Multi-axis slice.  Each argument corresponds to one axis and is either:
  //!   - an integer (fixes that axis, rank-reducing)
  //!   - All (keeps entire axis)
  //!   - Range (keeps sub-range along that axis)
  //! Example: t.slice(0, all, range(2, 5)) == numpy t[0, :, 2:5]
  template<typename First, typename... Rest>
  View<stored_type> slice(First first, Rest... rest)
  {
    auto strides = compute_strides();
    auto r = detail::compute_slice(shape_, strides, first, rest...);
    return {
      data_.data() + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }

  template<typename First, typename... Rest>
  View<const stored_type> slice(First first, Rest... rest) const
  {
    auto strides = compute_strides();
    auto r = detail::compute_slice(shape_, strides, first, rest...);
    return {
      data_.data() + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }

  //! Flat 1D view of all elements
  View<stored_type> flat()
  {
    return {data_.data(), {data_.size()}, {size_t(1)}};
  }
  View<const stored_type> flat() const
  {
    return {data_.data(), {data_.size()}, {size_t(1)}};
  }

  //--------------------------------------------------------------------------
  // Reductions and transforms

  //! Sum of all elements
  T sum() const
  {
    T s = T(0);
    for (size_t i = 0; i < data_.size(); ++i)
      s += data_[i];
    return s;
  }

  //! Sum along an axis, reducing rank by 1 (defined out-of-line below)
  Tensor<T> sum(size_t axis) const;

  //! Product of all elements
  T prod() const
  {
    T p = T(1);
    for (size_t i = 0; i < data_.size(); ++i)
      p *= data_[i];
    return p;
  }

  //! True if any element is nonzero
  bool any() const
  {
    for (size_t i = 0; i < data_.size(); ++i)
      if (data_[i])
        return true;
    return false;
  }

  //! True if all elements are nonzero
  bool all() const
  {
    for (size_t i = 0; i < data_.size(); ++i)
      if (!data_[i])
        return false;
    return true;
  }

  //! Flat index of the minimum element
  size_t argmin() const
  {
    return static_cast<size_t>(std::distance(data_.data(),
      std::min_element(data_.data(), data_.data() + data_.size())));
  }

  //! Reverse element order along an axis (e.g. flip(0) reverses rows)
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
  // Operators

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
  Tensor& operator+=(const Tensor& o)
  {
    for (size_t i = 0; i < data_.size(); ++i)
      data_[i] += o.data_[i];
    return *this;
  }
  Tensor& operator-=(const Tensor& o)
  {
    for (size_t i = 0; i < data_.size(); ++i)
      data_[i] -= o.data_[i];
    return *this;
  }
  Tensor& operator*=(const Tensor& o)
  {
    for (size_t i = 0; i < data_.size(); ++i)
      data_[i] *= o.data_[i];
    return *this;
  }
  Tensor& operator/=(const Tensor& o)
  {
    for (size_t i = 0; i < data_.size(); ++i)
      data_[i] /= o.data_[i];
    return *this;
  }

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
  Tensor operator*(const Tensor& o) const
  {
    Tensor r(shape_);
    for (size_t i = 0; i < data_.size(); ++i)
      r.data_[i] = data_[i] * o.data_[i];
    return r;
  }

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

  //! Compute row-major strides from shape
  vector<size_t> compute_strides() const
  {
    vector<size_t> strides(shape_.size());
    if (!shape_.empty()) {
      strides.back() = 1;
      for (int d = static_cast<int>(shape_.size()) - 2; d >= 0; --d)
        strides[d] = strides[d + 1] * shape_[d + 1];
    }
    return strides;
  }

  //--------------------------------------------------------------------------
  // Data members

  vector<size_t> shape_;
  vector<storage_type<T>> data_;
};

//==============================================================================
// Non-member operators (scalar op tensor)
//==============================================================================

template<typename T>
Tensor<T> operator*(T val, const Tensor<T>& arr)
{
  return arr * val;
}

template<typename T>
Tensor<T> operator+(T val, const Tensor<T>& arr)
{
  return arr + val;
}

// Mixed-type arithmetic: Tensor<T1> op Tensor<T2> -> Tensor<double>
// A SFINAE guard is used here, as without !is_same Tensor<T> * Tensor<T>
// would be ambiguous between the member operator* and this non-member function.
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

// Same SFINAE guard as operator* above.
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
// Out-of-line method definitions (require complete types)
//==============================================================================

template<typename T>
template<typename U>
View<T>& View<T>::operator=(const Tensor<U>& other)
{
  size_t n = size();
  for (size_t i = 0; i < n; ++i)
    data_[flat_to_offset(i)] = static_cast<T>(other.data()[i]);
  return *this;
}

template<typename T>
template<typename U>
View<T>& View<T>::operator+=(const Tensor<U>& o)
{
  size_t n = size();
  for (size_t i = 0; i < n; ++i)
    data_[flat_to_offset(i)] += o.data()[i];
  return *this;
}

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
// StaticTensor2D<T, R, C>: compile-time fixed 2D tensor.
//==============================================================================

template<typename T, size_t R, size_t C>
class StaticTensor2D {
public:
  using value_type = T;

  //--------------------------------------------------------------------------
  // Indexing

  //! Templated to accept enum class indices (e.g. GlobalTally, TallyResult)
  //! which don't implicitly convert to integer types.
  template<typename I0, typename I1>
  T& operator()(I0 i, I1 j)
  {
    return data_[static_cast<size_t>(i) * C + static_cast<size_t>(j)];
  }
  template<typename I0, typename I1>
  const T& operator()(I0 i, I1 j) const
  {
    return data_[static_cast<size_t>(i) * C + static_cast<size_t>(j)];
  }

  //--------------------------------------------------------------------------
  // Accessors

  T* data() { return data_; }
  const T* data() const { return data_; }
  constexpr size_t size() const { return R * C; }
  std::array<size_t, 2> shape() const { return {R, C}; }

  //--------------------------------------------------------------------------
  // Mutation

  void fill(T val) { std::fill(data_, data_ + R * C, val); }

  //--------------------------------------------------------------------------
  // Iterators

  T* begin() { return data_; }
  T* end() { return data_ + R * C; }
  const T* begin() const { return data_; }
  const T* end() const { return data_ + R * C; }

  //--------------------------------------------------------------------------
  // View accessors

  //! Multi-axis slice (same interface as Tensor/View).
  template<typename First, typename... Rest>
  View<T> slice(First first, Rest... rest)
  {
    vector<size_t> sh = {R, C};
    vector<size_t> st = {C, 1};
    auto r = detail::compute_slice(sh, st, first, rest...);
    return {data_ + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }
  template<typename First, typename... Rest>
  View<const T> slice(First first, Rest... rest) const
  {
    vector<size_t> sh = {R, C};
    vector<size_t> st = {C, 1};
    auto r = detail::compute_slice(sh, st, first, rest...);
    return {data_ + r.ptr_offset, std::move(r.shape), std::move(r.strides)};
  }

  //! Flat view (1D, contiguous)
  View<T> flat() { return {data_, {R * C}, {size_t(1)}}; }
  View<const T> flat() const { return {data_, {R * C}, {size_t(1)}}; }

private:
  //--------------------------------------------------------------------------
  // Data members

  T data_[R * C] = {};
};

//==============================================================================
// Non-member functions
//==============================================================================

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

template<typename T>
Tensor<T> ones(std::initializer_list<size_t> shape)
{
  vector<size_t> s(shape);
  return Tensor<T>(std::move(s), T(1));
}

template<typename T>
Tensor<T> ones(const vector<size_t>& shape)
{
  return Tensor<T>(shape, T(1));
}

template<typename T>
Tensor<T> zeros_like(const Tensor<T>& o)
{
  return Tensor<T>(o.shape(), T(0));
}

template<typename T, typename V>
Tensor<T> full_like(const Tensor<T>& o, V val)
{
  return Tensor<T>(o.shape(), static_cast<T>(val));
}

//! Return a 1D tensor of n evenly spaced values from start to stop (inclusive)
template<typename T>
Tensor<T> linspace(T start, T stop, size_t n)
{
  Tensor<T> result({n});
  if (n < 2) {
    result[0] = start;
    return result;
  }
  for (size_t i = 0; i < n; ++i) {
    result[i] =
      start + static_cast<T>(i) * (stop - start) / static_cast<T>(n - 1);
  }
  return result;
}

//! Concatenate two 1D tensors end-to-end
template<typename T>
Tensor<T> concatenate(const Tensor<T>& a, const Tensor<T>& b)
{
  size_t total = a.size() + b.size();
  Tensor<T> result({total});
  std::copy(a.data(), a.data() + a.size(), result.data());
  std::copy(b.data(), b.data() + b.size(), result.data() + a.size());
  return result;
}

//! Element-wise natural logarithm
template<typename T>
Tensor<T> log(const Tensor<T>& a)
{
  Tensor<T> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::log(a.data()[i]);
  return r;
}

//! Element-wise absolute value
template<typename T>
Tensor<T> abs(const Tensor<T>& a)
{
  Tensor<T> r(a.shape());
  for (size_t i = 0; i < a.size(); ++i)
    r.data()[i] = std::abs(a.data()[i]);
  return r;
}

//! Element-wise conditional: select from true_val where cond is true,
//! otherwise use false_val
template<typename T, typename V>
Tensor<T> where(
  const Tensor<bool>& cond, const Tensor<T>& true_val, V false_val)
{
  Tensor<T> r(cond.shape());
  for (size_t i = 0; i < cond.size(); ++i)
    r.data()[i] =
      cond.data()[i] ? true_val.data()[i] : static_cast<T>(false_val);
  return r;
}

//! Replace NaN/Inf values with finite substitutes
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

//! Type trait that is true for Tensor and StaticTensor2D.
//! Used by hdf5_interface.h to select the correct write_dataset overload.
template<typename T>
struct is_tensor : std::false_type {};

template<typename T>
struct is_tensor<Tensor<T>> : std::true_type {};

template<typename T, size_t R, size_t C>
struct is_tensor<StaticTensor2D<T, R, C>> : std::true_type {};

} // namespace tensor
} // namespace openmc

#endif // OPENMC_TENSOR_H
