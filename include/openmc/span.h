#ifndef OPENMC_SPAN_H
#define OPENMC_SPAN_H
#include <cstddef>   // for std::size_t, std::ptrdiff_t
#include <iterator>  // for std::begin, std::end
#include <stdexcept> // for std::out_of_range
#include <type_traits>

#include "openmc/vector.h"

namespace openmc {

template<typename T>
class span {
public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using iterator = T*;
  using const_iterator = const T*;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  /**
   * @brief Default constructor for an empty span.
   */
  span() noexcept : data_(nullptr), size_(0) {}

  /**
   * @brief Constructs a span from a pointer and size.
   * @param ptr Pointer to the first element.
   * @param count Number of elements in the span.
   */
  span(pointer ptr, size_type count) : data_(ptr), size_(count) {}

  /**
   * @brief Constructs a span from two pointers marking the span range.
   * @param first Pointer to the first element.
   * @param last Pointer past the last element.
   * @throws std::out_of_range if last < first.
   */
  span(pointer first, pointer last) : data_(first), size_(last - first)
  {
    if (last < first) {
      throw std::out_of_range("span: last pointer is before first pointer");
    }
  }

  /**
   * @brief Constructs a span from a non-const std::vector.
   * @param vec Reference to the vector to create a span from.
   */
  span(std::vector<T>& vec) : data_(vec.data()), size_(vec.size()) {}

  /**
   * @brief Constructs a span from a const std::vector.
   *
   * This is handling the semantics that a span<const double> is used
   * for read-only access into a vector.
   * @param vec Reference to the const vector to create a span from.
   */
  template<typename U = T,
    typename = std::enable_if_t<std::is_same<U, const T>::value>>
  span(const std::vector<std::remove_const_t<U>>& vec)
    : data_(vec.data()), size_(vec.size())
  {}

  /**
   * @brief Constructs a read-only span<const T> from a non-const span<T>.
   */
  template<typename U = T, typename = std::enable_if_t<std::is_const<U>::value>>
  span(const span<std::remove_const_t<U>>& other) noexcept
    : data_(other.data()), size_(other.size())
  {}

  /**
   * @brief Access an element without bounds checking.
   * @param index Index of the element to access.
   * @return Reference to the accessed element.
   */
  reference operator[](size_type index) { return data_[index]; }

  /**
   * @brief Access an element without bounds checking (const version).
   * @param index Index of the element to access.
   * @return Const reference to the accessed element.
   */
  const_reference operator[](size_type index) const { return data_[index]; }

  /**
   * @brief Access an element with bounds checking.
   * @param index Index of the element to access.
   * @return Reference to the accessed element.
   * @throws std::out_of_range if index is out of range.
   */
  reference at(size_type index)
  {
    if (index >= size_) {
      throw std::out_of_range("span: index out of range");
    }
    return data_[index];
  }

  /**
   * @brief Access an element with bounds checking (const version).
   * @param index Index of the element to access.
   * @return Const reference to the accessed element.
   * @throws std::out_of_range if index is out of range.
   */
  const_reference at(size_type index) const
  {
    if (index >= size_) {
      throw std::out_of_range("span: index out of range");
    }
    return data_[index];
  }

  /**
   * @brief Get a pointer to the underlying data.
   * @return Pointer to the data, or nullptr if the span is empty.
   */
  pointer data() noexcept { return data_; }

  /**
   * @brief Get a const pointer to the underlying data.
   * @return Const pointer to the data, or nullptr if the span is empty.
   */
  const_pointer data() const noexcept { return data_; }

  /**
   * @brief Get the number of elements in the span.
   * @return The size of the span.
   */
  size_type size() const noexcept { return size_; }

  /**
   * @brief Check if the span is empty.
   * @return True if the span is empty, false otherwise.
   */
  bool empty() const noexcept { return size_ == 0; }

  /**
   * @brief Get an iterator to the beginning of the span.
   * @return Iterator pointing to the first element.
   */
  iterator begin() noexcept { return data_; }

  /**
   * @brief Get a const iterator to the beginning of the span.
   * @return Const iterator pointing to the first element.
   */
  const_iterator begin() const noexcept { return data_; }

  /**
   * @brief Get a const iterator to the beginning of the span.
   * @return Const iterator pointing to the first element.
   */
  const_iterator cbegin() const noexcept { return data_; }

  /**
   * @brief Get an iterator to the end of the span.
   * @return Iterator pointing past the last element.
   */
  iterator end() noexcept { return data_ + size_; }

  /**
   * @brief Get a const iterator to the end of the span.
   * @return Const iterator pointing past the last element.
   */
  const_iterator end() const noexcept { return data_ + size_; }

  /**
   * @brief Get a const iterator to the end of the span.
   * @return Const iterator pointing past the last element.
   */
  const_iterator cend() const noexcept { return data_ + size_; }

  /**
   * @brief Access the first element.
   * @return Reference to the first element.
   * @throws std::out_of_range if the span is empty.
   */
  reference front()
  {
    if (empty()) {
      throw std::out_of_range("span::front(): span is empty");
    }
    return data_[0];
  }

  /**
   * @brief Access the first element (const version).
   * @return Const reference to the first element.
   * @throws std::out_of_range if the span is empty.
   */
  const_reference front() const
  {
    if (empty()) {
      throw std::out_of_range("span::front(): span is empty");
    }
    return data_[0];
  }

  /**
   * @brief Access the last element.
   * @return Reference to the last element.
   * @throws std::out_of_range if the span is empty.
   */
  reference back()
  {
    if (empty()) {
      throw std::out_of_range("span::back(): span is empty");
    }
    return data_[size_ - 1];
  }

  /**
   * @brief Access the last element (const version).
   * @return Const reference to the last element.
   * @throws std::out_of_range if the span is empty.
   */
  const_reference back() const
  {
    if (empty()) {
      throw std::out_of_range("span::back(): span is empty");
    }
    return data_[size_ - 1];
  }

private:
  pointer data_;
  size_type size_;
};

} // namespace openmc
#endif // OPENMC_SPAN_H
