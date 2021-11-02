#ifndef OPENMC_VECTOR_H
#define OPENMC_VECTOR_H

#include <algorithm> // for copy

namespace openmc {

template<typename T>
class vector {
public:
  // Types, aliases
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = T*;
  using const_iterator = const T*;


  // Constructors, destructors
  vector() : data_(nullptr), size_(0), capacity_(0) { }

  vector(size_type n) : vector() {
    this->reserve(n);
  }

  ~vector() {
    if (data_) {
#pragma omp target exit data map(delete: data_)
      delete [] data_;
    }
  }

  // Copy assignment
  vector& operator=(const vector& other);

  // Element access
  reference operator[](size_type pos) { return data_[pos]; }
  const_reference operator[](size_type pos) const { return data_[pos]; }
  reference front() { return data_[0]; }
  const_reference front() const { return data_[0]; }
  reference back() { return data_[size_ - 1]; }
  const_reference back() const { return data_[size_ - 1]; }
  pointer data() { return data_; }
  const_pointer data() const { return data_; }

  // Iterators
  iterator begin() { return data_; }
  const_iterator begin() const { return data_; }
  iterator end() { return data_ + size_; }
  const_iterator end() const { return data_ + size_; }
  const_iterator cbegin() const { return data_; }
  const_iterator cend() const { return data_ + size_; }

  // Capacity
  bool empty() const { return size_ == 0; }
  size_type size() const { return size_; }

  // Methods
  void clear() { size_ = 0; }

  void push_back(const T& value) {
    if (capacity_ == 0) {
      this->reserve(8);
    } else if (size_ == capacity_ ) {
      this->reserve(2*capacity_);
    }

    data_[size_] = value;
    ++size_;
  }


  void resize(size_type count) {
    this->reserve(count);
    size_ = count;
  }

  void reserve(size_type n) {
    if (n <= capacity_) return;

    // Make new allocation
    T* data_new = new T[n];

    // Copy existing elements
    for (size_type i = 0; i < size_; ++i) {
      data_new[i] = data_[i];
    }

    // Remove older allocation
    if (data_) {
      delete[] data_;
    }

    // Set data members
    data_ = data_new;
    capacity_ = n;
  }

  void copy_to_device() {
#pragma omp target enter data map(to: data_[:size_])
  }

  void copy_from_device() {
#pragma omp target exit data map(from: data_[:size_])
  }

  void release_device() {
#pragma omp target exit data map(release: data_[:size_])
  }

private:
  T* data_;
  size_type size_;
  size_type capacity_;
};

template<typename T>
vector<T>& vector<T>::operator=(const vector<T>& other)
{
  this->resize(other.size());
  std::copy(other.cbegin(), other.cend(), this->begin());
  return *this;
}

} // namespace openmc

#endif
