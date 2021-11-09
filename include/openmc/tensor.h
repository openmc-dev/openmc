#ifndef TENSOR_H
#define TENSOR_H

/*
 * In an implementation of OpenMC that offloads computations to an accelerator,
 * we may need to provide replacements for standard library containers and
 * algorithms that have no native implementations on the device of interest.
 * Because some developers are currently in the process of creating such code,
 * introducing the below typedef lessens the amount of rebase conflicts that
 * happen as they rebase their code on OpenMC's develop branch.
 *
 */

#include <iostream>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>

#include "openmc/error.h"
#include "openmc/memory.h"

#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define HOST __host__
#else
#define HOSTDEVICE
#define HOST
#endif

namespace openmc {

namespace impl {
// Use with enable_if to check if template argument makes sense
// as a tensor size definition. It has to be a vector-like container
// of integral values. Would be a lot nicer-looking with C++20
// concepts.

// C++17 provides an implementation of this
template<class...>
using void_t = void;

template<typename Vector, typename Filler = void>
struct describes_tensor_shape : std::false_type {};
template<typename Vector>
struct describes_tensor_shape<Vector,
  void_t<decltype(std::declval<Vector&>().size()),
    decltype(std::declval<Vector&>().data()),
    decltype(std::declval<Vector&>()[0]),
    typename std::enable_if<
      std::is_integral<typename Vector::value_type>::value>::type>>
  : std::true_type {};

// This creates a tuple of pre-specified length of a single type repeated N
// times.
using size_type = unsigned int;
template<typename Base, typename Seq>
struct expander;
template<typename Base, std::size_t... Is>
struct expander<Base, std::index_sequence<Is...>> {
  template<typename E, std::size_t>
  using elem = E;
  using type = std::tuple<elem<Base, Is>...>;
};

// bracket_indx is a functor that recursively generates functors at compile time
// in order to compute the tensor index without any loops for arbitrary tensor
// rank
template<int Idx, int Rank>
struct bracket_indx {
  using tuple_type =
    typename expander<size_type, std::make_index_sequence<Rank>>::type;
  HOSTDEVICE static size_type calc(
    size_type coeff, tuple_type indices, const size_type* size)
  {
    return coeff * std::get<Idx>(indices) +
           bracket_indx<Idx - 1, Rank>::calc(coeff * size[Idx], indices, size);
  }
};
template<int Rank>
struct bracket_indx<0, Rank> { // Specialize on base case to not recurse
  using tuple_type =
    typename expander<size_type, std::make_index_sequence<Rank>>::type;
  HOSTDEVICE static size_type calc(
    size_type coeff, tuple_type indices, const size_type* size)
  {
    return coeff * std::get<0>(indices);
  }
};
} // namespace impl

template<typename T, int Rank, typename Alloc = UnifiedAllocator<T>>
class tensor {
public:
  using value_type = T;
  using size_type = impl::size_type;
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
  size_type size_[Rank];
  size_type capacity_;
  Alloc alloc_;

public:
  HOSTDEVICE size_type num_elements() const
  {
    size_type result = 1;
#pragma unroll
    for (int i = 0; i < Rank; ++i) {
      result *= size_[i];
    }
    return result;
  }

  HOST tensor() : begin_(nullptr), capacity_(0)
  {
    for (int i = 0; i < Rank; ++i)
      size_[i] = 0;
  }

  HOST tensor(const size_type* new_size, const T& value = T())
  {
    reserve(new_size, value);
  }

  // Construction off of anything vector-like defining a tensor shape
  template<typename Vector,
    typename = std::enable_if_t<impl::describes_tensor_shape<Vector>::value>>
  HOST tensor(Vector const& new_size, const T& value = T())
  {
    if (new_size.size() != Rank)
      throw std::runtime_error("Tensor rank doesn't match input shape size.");
    reserve(new_size.data(), value);
  }

  // The copy constructor should create a totally independent deep copy
  HOST tensor(tensor const& copy_from)
  {
    begin_ = alloc_.allocate(copy_from.num_elements());
    for (int i = 0; i < Rank; ++i)
      size_[i] = copy_from.size(i);
    capacity_ = total_size(size_);
    alloc_ = copy_from.alloc_;
    for (size_type i = 0; i < capacity_; ++i)
      new (begin_ + i) T(copy_from.begin_[i]);
  }

  HOST tensor& operator=(tensor&& move_from)
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

  // Allow construction from tensor
  HOST tensor& operator=(tensor<T, Rank, Alloc> const& copy_from)
  {
    begin_ = alloc_.allocate(copy_from.num_elements());
    for (int i = 0; i < Rank; ++i)
      size_[i] = copy_from.size_[i];
    capacity_ = num_elements();
    alloc_ = copy_from.alloc_;
    for (size_type i = 0; i < size_; ++i)
      new (begin_ + i) T(copy_from[i]);
    return *this;
  }

  // Allow copy assignment from anything that looks like an xtensor
  template<typename View>
  HOST tensor& operator=(View const& view)
  {
    if (view.shape().size() != Rank)
      throw std::invalid_argument {
        "operator= template on openmc::tensor requires ranks to match!"};

    // NOTE: this is not guaranteed to work, since shape()
    // may not be contiguous storage.
    resize(&view.shape()[0]);
    auto num_el = num_elements();
    for (unsigned i = 0; i < num_el; ++i)
      begin_[i] = view(i);
    return *this;
  }

  // Allow division and multiplication of the array by scalars in the
  // way that you would expect. This allows drop-in replacement of xtensors
  // to a limited extent. Template expressions would make this
  // far better, but we actually don't use those much in OpenMC.
  HOST tensor& operator*=(double multiplier)
  {
    for (unsigned i = 0; i < size_; ++i)
      begin_[i] *= multiplier;
    return *this;
  }
  HOST tensor& operator/=(double divisor)
  {
    double multiplier = 1.0 / divisor;
    for (unsigned i = 0; i < size_; ++i)
      begin_[i] *= multiplier;
    return *this;
  }

  // The move constructor may need to run on the device in the case of
  // construction of polymorphic objects living on GPU that contain tensors.
#pragma hd_warning_disable
  HOSTDEVICE tensor(tensor&& move_from)
    : begin_(std::move(move_from.begin_)), size_(std::move(move_from.size_)),
      capacity_(std::move(move_from.capacity_)),
      alloc_(std::move(move_from.alloc_))
  {
    // Remove responsibility of the other one for the memory it held
    move_from.begin_ = nullptr;
    for (int i = 0; i < Rank; ++i)
      move_from.size_[i] = 0;
    move_from.capacity_ = 0;
  }

  HOST ~tensor()
  {
    auto num_elements = total_size(&size_[0]);
    for (size_type i = 0; i < num_elements; ++i)
      (begin_ + i)->~T();
    if (capacity_)
      alloc_.deallocate(begin_, capacity_);
  }

  // Resize off of anything vector-like defining a tensor shape
  template<typename Vector,
    typename = std::enable_if_t<impl::describes_tensor_shape<Vector>::value>>
  HOST void resize(Vector const& new_size, const T& default_value = T())
  {
    if (new_size.size() != Rank)
      fatal_error("Mismatch in tensor::resize and new size array");
    size_type num_elements = total_size(&new_size[0]);
    if (num_elements > capacity_)
      reserve(&new_size[0]);
    // set new things to be default_value
    for (size_type i = 0; i < num_elements; ++i)
      begin_[i] = default_value;
    for (int i = 0; i < Rank; ++i)
      size_[i] = new_size[i];
  }

  HOST void resize(const size_type* new_size, T const& default_value = T())
  {
    size_type num_elements = total_size(new_size);
    if (num_elements > capacity_)
      reserve(new_size);
    // set new things to be default_value
    for (size_type i = 0; i < num_elements; ++i)
      begin_[i] = default_value;
    for (int i = 0; i < Rank; ++i)
      size_[i] = new_size[i];
  }

  HOSTDEVICE bool empty() const { return num_elements() == 0; }
  HOSTDEVICE pointer data() const { return begin_; }
  HOSTDEVICE pointer data() { return begin_; }
  HOSTDEVICE size_type size(int i) const { return size_[i]; }
  template<int Idx>
  HOSTDEVICE size_type size() const
  {
    return size_[Idx];
  }
  HOSTDEVICE std::array<size_type, Rank> shape() const
  {
    std::array<size_type, Rank> result;
    for (int i = 0; i < Rank; ++i)
      result[i] = size_[i];
    return result;
  }
  HOSTDEVICE size_type capacity() const { return capacity_; }
  HOSTDEVICE iterator begin() { return begin_; }
  HOSTDEVICE iterator end() { return begin_ + num_elements(); }
  HOSTDEVICE iterator begin() const { return begin_; }
  HOSTDEVICE iterator end() const { return begin_ + num_elements(); }
  HOSTDEVICE const_iterator cbegin() const { return begin_; }
  HOSTDEVICE const_iterator cend() const { return begin_ + num_elements(); }
  HOST reverse_iterator rbegin() const
  {
    return std::make_reverse_iterator(begin_ + num_elements());
  }
  HOST reverse_iterator rend() const
  {
    return std::make_reverse_iterator(begin_);
  }
  HOST const_reverse_iterator crbegin() const
  {
    return std::make_reverse_iterator(begin_ + num_elements());
  }
  HOST const_reverse_iterator crend() const
  {
    return std::make_reverse_iterator(begin_);
  }
  HOSTDEVICE const_reference back() const { return begin_[num_elements() - 1]; }
  HOSTDEVICE const_reference front() const { return begin_[0]; }
  HOSTDEVICE reference back() { return begin_[num_elements() - 1]; }
  HOSTDEVICE reference front() { return begin_[0]; }

  // Enable host-device synchronization functions if replicated memory is being
  // used
  // template<typename U = pointer>
  // std::enable_if_t<std::is_same<U, DualPointer<T>>::value> syncToDevice()
  // {
  //   cudaMemcpy(begin_.device_pointer(), begin_.host_pointer(),
  //     sizeof(value_type) * size_, cudaMemcpyHostToDevice);
  // }
  // template<typename U = pointer>
  // std::enable_if_t<std::is_same<U, DualPointer<T>>::value> syncToHost()
  // {
  //   cudaMemcpy(begin_.host_pointer(), begin_.device_pointer(),
  //     sizeof(value_type) * size_, cudaMemcpyDeviceToHost);
  // }

  template<typename... Args>
  HOSTDEVICE T& operator()(Args... args)
  {
    constexpr size_type one = 1;
    return begin_[impl::bracket_indx<Rank - 1, Rank>::calc(
      one, std::forward_as_tuple(std::forward<Args>(args)...), size_)];
  }
  template<typename... Args>
  HOSTDEVICE T& operator()(Args... args) const
  {
    constexpr size_type one = 1;
    return begin_[impl::bracket_indx<Rank - 1, Rank>::calc(
      one, std::forward_as_tuple(std::forward<Args>(args)...), size_)];
  }

private:
  HOSTDEVICE size_type total_size(const size_type* values)
  {
    size_type result = 1;
    for (int i = 0; i < Rank; ++i) {
      result *= values[i];
    }
    return result;
  }

  // Construct, default-initializing elements in a tensor of
  // length new_size. Does no checking on the argument.
  HOST void reserve(const size_type* new_size, const T& value = T())
  {
    size_type num_elements = 1;
    for (int i = 0; i < Rank; ++i)
      num_elements *= new_size[i];

    begin_ = alloc_.allocate(num_elements);
    capacity_ = num_elements;
    for (size_type i = 0; i < num_elements; ++i)
      new (begin_ + i) T(value);
    for (int i = 0; i < Rank; ++i)
      size_[i] = new_size[i];
  }
};

template<typename T, int Rank, typename Alloc>
HOSTDEVICE bool operator==(
  const tensor<T, Rank, Alloc>& first, const tensor<T, Rank, Alloc>& second)
{
  if (first.size() != second.size())
    return false;
  for (typename tensor<T, Rank, Alloc>::size_type i = 0; i < first.size(); ++i)
    if (first[i] != second[i])
      return false;
  return true;
}

template<typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, tensor<T, 1, Alloc> const& output)
{
  auto num_elem = output.size(0);
  using size_type = typename tensor<T, 1, Alloc>::size_type;
  for (size_type i = 0; i < num_elem; ++i)
    os << output(i) << std::endl;
  return os;
}
template<typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, tensor<T, 2, Alloc> const& output)
{
  auto num_row = output.size(0);
  auto num_col = output.size(1);
  using size_type = typename tensor<T, 2, Alloc>::size_type;
  for (size_type i = 0; i < num_row; ++i) {
    for (size_type j = 0; j < num_col - 1; ++j)
      os << output(i, j) << " ";
    os << output(i, num_col - 1) << std::endl;
  }
  return os;
}

} // namespace openmc
#endif // TENSOR_H
