#ifndef OPENMC_MEMORY_H
#define OPENMC_MEMORY_H

#include <cstddef>
#include <memory>
#include <type_traits>

#ifdef __CUDACC__
#include <cub/cub.cuh>
#endif

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
    cudaError_t error_code = cudaMallocManaged(&tmp, n * sizeof(value_type));
    if (error_code != cudaSuccess)
      CubDebugExit(error_code);
    return tmp;
  }
  inline void deallocate(pointer p, size_type)
  {
    if (p == nullptr) return;
    cudaError_t error_code = cudaFree(p);
    // For the latter, this happens when running deallocate in destructors
    // in on of OpenMC's many global variables which destruct after main().
    // The CUDA RT API has already handled freeing the memory in that case,
    // so that error is safe to ignore.
    if (error_code != cudaSuccess && error_code != cudaErrorCudartUnloading)
      CubDebugExit(error_code);
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

// Expressive way to give a yes/no answer regarding whether
// pinned memory is to be used on the host.
enum class PinnedAllocationOnHost { yes, no };

template<typename T, PinnedAllocationOnHost UsePinned>
class ReplicatedAllocator;

// Keeps track of replicated memory, and acts like a normal pointer as
// appropriate
template<typename T>
class DualPointer {
private:
  T* host;
  T* device;

public:
  __host__ __device__ T& operator*()
  {
#ifdef __CUDA_ARCH__
    return *device;
#else
    return *host;
#endif
  }
  __host__ __device__ T const& operator*() const
  {
#ifdef __CUDA_ARCH__
    return *device;
#else
    return *host;
#endif
  }
  __host__ __device__ operator T*()
  {
#ifdef __CUDA_ARCH__
    return device;
#else
    return host;
#endif
  }
  __host__ __device__ operator T const *() const
  {
#ifdef __CUDA_ARCH__
    return device;
#else
    return host;
#endif
  }

  __host__ __device__ T const& operator[](unsigned int indx) const
  {
#ifdef __CUDA_ARCH__
    return device[indx];
#else
    return host[indx];
#endif
  }
  __host__ __device__ T& operator[](unsigned int indx)
  {
#ifdef __CUDA_ARCH__
    return device[indx];
#else
    return host[indx];
#endif
  }

  __host__ __device__ T* operator+(unsigned int offset)
  {
#ifdef __CUDA_ARCH__
    return device + offset;
#else
    return host + offset;
#endif
  }
  __host__ __device__ T const* operator+(unsigned int offset) const
  {
#ifdef __CUDA_ARCH__
    return device + offset;
#else
    return host + offset;
#endif
  }

  __host__ __device__ DualPointer() : host(nullptr), device(nullptr) {}

  // Copy a given pointer into each. This can be used with managed memory
  __host__ __device__ DualPointer(T* const& ptr) : host(ptr), device(ptr) {}
  __host__ __device__ DualPointer& operator=(T* const& ptr)
  {
    host = ptr;
    device = ptr;
    return *this;
  }

  // These should exclusively be used for cudaMemcpy calls
  __host__ T* host_pointer() { return host; }
  __host__ T* device_pointer() { return device; }

  DualPointer(DualPointer const&) = default;
  DualPointer& operator=(DualPointer const&) = default;
  DualPointer(DualPointer&&) = default;
  DualPointer& operator=(DualPointer&&) = default;
  ~DualPointer() = default;

  friend class ReplicatedAllocator<T, PinnedAllocationOnHost::no>;
  friend class ReplicatedAllocator<T, PinnedAllocationOnHost::yes>;
};

template<typename T,
  PinnedAllocationOnHost UsePinned = PinnedAllocationOnHost::no>
class ReplicatedAllocator {
public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef DualPointer<T> pointer;

  // Not const-correct here. Not sure how to easily implement that in
  // DualPointer.
  typedef DualPointer<T> const_pointer;

  typedef T& reference;
  typedef const T& const_reference;

  inline ReplicatedAllocator() throw() {}
  template<typename T2>
  inline ReplicatedAllocator(const ReplicatedAllocator<T2>&) throw()
  {}
  inline ~ReplicatedAllocator() throw() {}
  inline pointer address(reference r) { return &r; }
  inline const_pointer address(const_reference r) const { return &r; }
  inline pointer allocate(size_type n)
  {
    pointer tmp;

    if (UsePinned == PinnedAllocationOnHost::yes) {
      cudaError_t error_code =
        cudaMallocHost(&tmp.host, n * sizeof(value_type));
      if (error_code != cudaSuccess)
        CubDebugExit(error_code);
    } else {
      tmp.host = (T*)malloc(n * sizeof(value_type));
    }

    // Allocate on device
    cudaError_t error_code = cudaMalloc(&tmp.device, n * sizeof(value_type));
    if (error_code != cudaSuccess)
      CubDebugExit(error_code);

    return tmp;
  }
  inline void deallocate(pointer p, size_type)
  {
    if (p.host) {
      if (UsePinned == PinnedAllocationOnHost::yes) {
        cudaError_t error_code = cudaFree(p.host);
        if (error_code != cudaSuccess && error_code != cudaErrorCudartUnloading)
          CubDebugExit(error_code);
      } else {
        free(p.host);
      }
    }

    // See comment in UnifiedAllocator on this method
    if (p.device) {
      cudaError_t error_code = cudaFree(p.device);
      if (error_code != cudaSuccess && error_code != cudaErrorCudartUnloading)
        CubDebugExit(error_code);
    }
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
    typedef ReplicatedAllocator<T2> other;
  };
  bool operator!=(const ReplicatedAllocator<T>& other) const
  {
    return !(*this == other);
  }
  bool operator==(const ReplicatedAllocator<T>& /* other */) const
  {
    return true;
  }
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
__global__ void run_move_constructor_on_device(T* dist)
{
  static_assert(std::is_move_constructible<T>::value,
    "Polymorphic objects to be put on device must be move constructible.");

  // NOTE: possibly need to use intermediate temporary object here.
  // Seems to work though. I think this is not necessarily guaranteed
  // to work though.
  new (dist) T(std::move(*dist));
}

// Allocates unified memory if running in GPU mode, else acts
// like operator new.
template<typename T, typename... Args>
T* unified_new(Args... args)
{
  T* loc = nullptr;

  cudaMallocManaged(&loc, sizeof(T));

  // Run the constructor as usual on the host
  new (loc) T(args...);
  return loc;
}

// Returns a pointer to the same object constructed on both the host,
// and the device separately. This is necessary for polymorphic stuff.
template<typename T, typename... Args>
std::pair<T*, T*> replicated_new(Args... args)
{
  T* loc_host = nullptr;
  T* loc_device = nullptr;

  loc_host = static_cast<T*>(malloc(sizeof(T)));
  cudaMalloc(&loc_device, sizeof(T));

  // Run the constructor as usual on the host
  new (loc_host) T(args...);

  // Now, copy the host data to device, and run the copy constructor there.
  // This will correctly initialize vtables on the device. Without re-running
  // a constructor on the device, polymorphism fails.
  cudaMemcpy(loc_device, loc_host, sizeof(T), cudaMemcpyHostToDevice);
  run_move_constructor_on_device<<<1, 1>>>(loc_device);

  return {loc_host, loc_device};
}

/**
 * This is a partial implementation of std::unique_ptr that can be
 * used on the GPU. By default, CUDA does not define __device__ functions
 * for this, even though they are just as efficient as raw pointers on
 * the device. This fills that gap. Memory is assumed to be managed by the
 * host and made available on the device by use of the unified memory
 * capability.
 *
 * TODO: this could be cleaned up to internally use the DualPointer
 * class rather than two raw pointers.
 */
template<class T>
class unique_ptr {
public:
  typedef T* pointer;
  typedef T element_type;

  // This one assumes that the pointer was allocated
  // using cudaMallocManaged, and thus contains a pointer usable
  // from both host and device.
  explicit unique_ptr(pointer p) : ptr(p), ptr_dev(p) {}

  // This one assumes the pointers were separately allocated
  // and are replications of the same data. This is done so that
  // both copies of host and device compatible vtables are kept
  // for polymorphic data. make_replicated handles this.
  explicit unique_ptr(std::pair<pointer, pointer> ptr_pair)
    : ptr(ptr_pair.first), ptr_dev(ptr_pair.second)
  {}

  unique_ptr(unique_ptr&& u) { std::tie(ptr, ptr_dev) = u.release(); }

  template<class U>
  unique_ptr(unique_ptr<U>&& u)
  {
    std::tie(ptr, ptr_dev) = u.release();
  }

  // Constructors/operators that unique_ptr must rid of to attain its desired
  // functionality
  unique_ptr(unique_ptr const&) = delete;
  unique_ptr& operator=(unique_ptr const&) = delete;
  unique_ptr& operator=(unique_ptr&& u)
  {
    free_mem(); // free memory currently held, if any
    std::tie(ptr, ptr_dev) = u.release();
    return *this;
  }

  void free_mem()
  {
    // If NOT using managed memory, must separately free device memory
    if (ptr_dev != ptr) {
      cudaFree(ptr_dev);
    }
    // Now, free whatever lives on the host OR the managed memory
    // if that was not called in the prior cudaFree.
    if (ptr != nullptr) {
      ptr->~T();
      cudaFree(ptr);
    }
  }

  ~unique_ptr() { free_mem(); }

  // observers
  __host__ __device__ typename std::add_lvalue_reference<T>::type operator*()
    const
  {
#ifdef __CUDA_ARCH__
    return *ptr_dev;
#else
    return *ptr;
#endif
  }

  // observers
  __host__ __device__ typename std::add_lvalue_reference<T>::type operator*()
  {
#ifdef __CUDA_ARCH__
    return *ptr_dev;
#else
    return *ptr;
#endif
  }

  // These observers can be called on the GPU or the CPU.
  __host__ __device__ pointer operator->() const noexcept
  {
#ifdef __CUDA_ARCH__
    return ptr_dev;
#else
    return ptr;
#endif
  }

  // These observers can be called on the GPU or the CPU.
  __host__ __device__ pointer operator->() noexcept
  {
#ifdef __CUDA_ARCH__
    return ptr_dev;
#else
    return ptr;
#endif
  }

  __host__ __device__ pointer get() const noexcept
  {
#ifdef __CUDA_ARCH__
    return ptr_dev;
#else
    return ptr;
#endif
  }

  __host__ __device__ pointer get_dev() const noexcept { return ptr_dev; }
  __host__ __device__ pointer get_host() const noexcept { return ptr; }

  std::tuple<pointer, pointer> release() noexcept
  {
    pointer tmp = ptr;
    pointer tmp_dev = ptr_dev;
    ptr = nullptr;
    ptr_dev = nullptr;
    return {tmp, tmp_dev};
  }

  __host__ __device__ explicit operator bool() const noexcept
  {
#ifdef __CUDA_ARCH__
    return ptr_dev != nullptr;
#else
    return ptr != nullptr;
#endif
  }

  // The following two methods synchronize data between host and device manually
  // if a polymorphic object is being pointed to by this unique_ptr object.
  // Because polymorphic things can't use managed memory (the vtable would not
  // be set up for each host and device), these are done manually for now. This
  // could be made automatic with some ease.
  __host__ void syncHostToDevice()
  {
    if (!ptr_dev)
      return; // don't do anything if not using managed memory
    cudaError_t error_code =
      cudaMemcpy(ptr_dev, ptr, sizeof(T), cudaMemcpyHostToDevice);
    if (error_code != cudaSuccess)
      CubDebugExit(error_code);
  }
  __host__ void syncDeviceToHost()
  {
    if (!ptr_dev)
      return; // don't do anything if not using managed memory
    cudaError_t error_code =
      cudaMemcpy(ptr, ptr_dev, sizeof(T), cudaMemcpyDeviceToHost);
    if (error_code != cudaSuccess)
      CubDebugExit(error_code);
  }

private:
  pointer ptr {nullptr};
  pointer ptr_dev {nullptr};
};

/**
 * This is std::make_unique, but it does not use the "new" operator.
 * Rather, it will call a function which optionally uses CUDA
 * unified memory if we're running in GPU mode. SFINAE selects whether
 * the replicated or managed memory approach is taken depending on
 * whether the object under consideration is polymorphic.
 */
template<typename T, typename... Args>
std::enable_if_t<std::is_polymorphic<T>::value, unique_ptr<T>> make_unique(
  Args&&... args)
{
  return unique_ptr<T>(replicated_new<T>(std::forward<Args>(args)...));
}
template<typename T, typename... Args>
std::enable_if_t<!std::is_polymorphic<T>::value, unique_ptr<T>> make_unique(
  Args&&... args)
{
  return unique_ptr<T>(unified_new<T>(std::forward<Args>(args)...));
}

// Below are the expected standard library routines and data if not compiling
// for GPU.

#else

// template<class T, class D = std::default_delete<T>>
// using unique_ptr = std::unique_ptr<T, D>;
using std::make_unique;
using std::unique_ptr;

#endif

} // namespace openmc

#endif // OPENMC_MEMORY_H
