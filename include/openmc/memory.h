#ifndef OPENMC_MEMORY_H
#define OPENMC_MEMORY_H

#include <cstddef>
#include <memory>
#include <type_traits>

#ifdef __CUDACC__
#include <cub/cub.cuh>
#endif

#include "openmc/settings.h"

namespace openmc {

#ifdef __CUDACC__

/**
 * Checks if OpenMC is being run in GPU mode. If it isn't, this
 * allocats memory like a usual allocator. If it is, it allocates
 * CUDA unified memory which is virtually available on both host
 * and the device.
 */
template<typename T>
class UnifiedAllocator {
public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;

  inline UnifiedAllocator() throw() {}
  template<typename T2>
  inline UnifiedAllocator(const UnifiedAllocator<T2>&) throw()
  {}
  inline ~UnifiedAllocator() throw() {}
  inline pointer address(reference r) { return &r; }
  inline const_pointer address(const_reference r) const { return &r; }
  inline pointer allocate(size_type n)
  {
    pointer tmp;
    if (settings::gpu_mode) {
      cudaError_t error_code =
        cudaMallocManaged(&tmp, n * sizeof(value_type), cudaMemAttachGlobal);
      if (error_code != cudaSuccess)
        CubDebugExit(error_code);
    } else
      tmp = (pointer)malloc(n * sizeof(value_type));
    return tmp;
  }
  inline void deallocate(pointer p, size_type)
  {
    if (settings::gpu_mode) {
      cudaError_t error_code = cudaFree(p);
      if (error_code != cudaSuccess)
        CubDebugExit(error_code);
    } else
      free(p);
  }

  // Need to be able to construct in the newly allocated memory in a few ways
  inline void construct(pointer p) { new (p) value_type(); }
  inline void construct(pointer p, const value_type& wert)
  {
    new (p) value_type(wert);
  }
  inline void construct(pointer p, value_type& wert)
  {
    new (p) value_type(wert);
  }

  // This method is required for constructing vectors of unique_ptrs
  template<class U, class... Args>
  inline void construct(U* p, Args&&... args)
  {
    new (p) U(std::forward<Args>(args)...);
  }

  inline void destroy(pointer p) { p->~value_type(); }
  inline size_type max_size() const throw()
  {
    return size_type(-1) / sizeof(value_type);
  }
  template<typename T2>
  struct rebind {
    typedef UnifiedAllocator<T2> other;
  };
  bool operator!=(const UnifiedAllocator<T>& other) const
  {
    return !(*this == other);
  }
  bool operator==(const UnifiedAllocator<T>& /* other */) const { return true; }
};

// This can be used for running constructors on the device using arguments from
// the host. You have to note that stuff like XML parsing for the constructor
// won't work, in this case.
template<typename T, typename... Args>
__global__ void run_constructor_on_device(T* dist, Args... args)
{
  new (dist) T(args...);
}

template<typename T>
__global__ void run_copy_constructor_on_device(T* dist)
{
  static_assert(std::is_copy_constructible<T>::value,
    "Polymorphic objects to be used on the GPU must be copy constructible.");
  new (dist) T(*dist);
}

// Allocates unified memory if running in GPU mode, else acts
// like operator new.
template<typename T, typename... Args>
T* openmc_new(Args... args)
{
  T* loc = nullptr;

  if (settings::gpu_mode)
    cudaMallocManaged(&loc, sizeof(T));
  else
    loc = (T*)malloc(sizeof(T));

  // Run the constructor as usual on the host
  new (loc) T(args...);

  // Now, copy the host data to device, and run the copy constructor there.
  // This will correctly initialize vtables on the device. Without re-running
  // a constructor on the device, polymorphism fails.
  if (std::is_polymorphic<T>::value && settings::gpu_mode)
    run_copy_constructor_on_device<<<1, 1>>>(loc);

  return loc;
}

/**
 * This is a partial implementation of std::unique_ptr that can be
 * used on the GPU. By default, CUDA does not define __device__ functions
 * for this, even though they are just as efficient as raw pointers on
 * the device. This fills that gap. Memory is assumed to be managed by the
 * host and made available on the device by use of the unified memory
 * capability.
 */
template<class T>
class unique_ptr {
public:
  typedef T* pointer;
  typedef T element_type;

  // Allowable constructors
  explicit unique_ptr(pointer p) : ptr(p) {}

  unique_ptr(unique_ptr&& u) : ptr(u.release()) {}

  template<class U>
  unique_ptr(unique_ptr<U>&& u) : ptr(u.release())
  {}

  // Constructors/operators that unique_ptr must rid of to attain its desired
  // functionality
  unique_ptr(unique_ptr const&) = delete;
  unique_ptr& operator=(unique_ptr const&) = delete;

  ~unique_ptr()
  {
    if (ptr != nullptr)
      cudaFree(ptr);
  }

  // observers
  __host__ __device__ typename std::add_lvalue_reference<T>::type operator*()
    const
  {
    return *ptr;
  }

  // These observers can possibly be called on the GPU
  __host__ __device__ pointer operator->() const noexcept { return ptr; }

  __host__ __device__ pointer get() const noexcept { return ptr; }

  pointer release() noexcept
  {
    pointer tmp = ptr;
    ptr = nullptr;
    return tmp;
  }

  __host__ __device__ explicit operator bool() const noexcept
  {
    return ptr != nullptr;
  }

private:
  pointer ptr {nullptr};
};

/**
 * This is std::make_unique, but it does not use the "new" operator.
 * Rather, it will call a function which optionally uses CUDA
 * unified memory if we're running in GPU mode.
 */
template<typename T, typename... Args>
unique_ptr<T> make_unique(Args&&... args)
{
  return unique_ptr<T>(openmc_new<T>(std::forward<Args>(args)...));
}

// Below are the expected standard library routines and data if not compiling
// for GPU.

#else

template<class T, class D>
using unique_ptr = std::unique_ptr<T, D>;

template<class T, class... Args>
using make_unique = std::make_unique<T, Args...>

#endif

} // namespace openmc

#endif // OPENMC_MEMORY_H
