#ifdef SYCL_SORT
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <sycl/sycl.hpp>
#include <omp.h>
#include <map>
#include <iostream>

// All code in this file (except for the last function) is boilerplate and
// was provided by Thomas Applencourt.

namespace sycl = cl::sycl;

struct syclDeviceInfo {
  sycl::context sycl_context;
  sycl::device sycl_device;
  ze_context_handle_t ze_context;
};

// General case where each Context can have target multiple Device
std::vector<struct syclDeviceInfo> xomp_get_infos_devices() {
  static std::vector<struct syclDeviceInfo> ompDeviceId2Context;
  if (!ompDeviceId2Context.empty())
    return ompDeviceId2Context;

  // This datascructure is a vector used to map a OpenMP device ID to a sycl context and device
  ompDeviceId2Context.resize(omp_get_num_devices());
  //   Map each level zero (pcontext) to a vector a sycl::device.
  //   This is requied to create SYCL::context spaming multiple devices.
  std::map<ze_context_handle_t, std::vector<sycl::device>> hContext2device;
  for (int D=0; D< omp_get_num_devices(); D++) {
    omp_interop_t o = 0;
    #pragma omp interop init(targetsync: o) device(D)
    int err = -1;
    const ze_driver_handle_t hPlatform = static_cast<ze_driver_handle_t>(omp_get_interop_ptr(o, omp_ipr_platform, &err));
    assert (err >= 0 && "omp_get_interop_ptr(omp_ipr_platform)");
    const ze_context_handle_t hContext = static_cast<ze_context_handle_t>(omp_get_interop_ptr(o, omp_ipr_device_context, &err));
    assert (err >= 0 && "omp_get_interop_ptr(omp_ipr_device_context)");
    const ze_device_handle_t hDevice =  static_cast<ze_device_handle_t>(omp_get_interop_ptr(o, omp_ipr_device, &err));
    assert (err >= 0 && "omp_get_interop_ptr(omp_ipr_device)");
    #pragma omp interop destroy(o)

#ifdef MAKE_PLATFORM
    const sycl::platform sycl_platform = sycl::make_platform<sycl::backend::ext_oneapi_level_zero>(hPlatform);
#else
    const sycl::platform sycl_platform = sycl::ext::oneapi::level_zero::make<sycl::platform>(hPlatform);
#endif

#ifdef MAKE_DEVICE
    const sycl::device sycl_device = sycl::make_device<sycl::backend::ext_oneapi_level_zero>(hDevice);
#else
    const sycl::device sycl_device = sycl::level_zero::make<sycl::device>(sycl_platform, hDevice);
#endif
    // Store the Level_zero context assciated to this OMP Device ID. This will be required to create the SYCL context latter
    ompDeviceId2Context[D].ze_context = hContext;
    // Store the sycl device associate to this  OMP Device ID.
    ompDeviceId2Context[D].sycl_device = sycl_device;

    // Now we need to associate the level zero Context to the  device
    hContext2device[hContext].push_back(sycl_device);
  }

  // Construct sycl::contexts who stawn multiple openmp device, if possible.
  // This is N2, but trivial to make it log(N)
  for ( const auto& [ hContext, sycl_devices]: hContext2device ) {
#ifdef MAKE_CONTEXT
    // One more layer of fun (https://github.com/intel/llvm-test-suite/blob/intel/SYCL/Plugin/interop-level-zero.cpp)
    // No more "KeepOwnership" for context so the code segfault...
    sycl::backend_input_t<sycl::backend::ext_oneapi_level_zero, sycl::context> hContextInteropInput = {hContext, sycl_devices};
    const sycl::context sycl_context = sycl::make_context<sycl::backend::ext_oneapi_level_zero>(hContextInteropInput);
#else
    const sycl::context sycl_context = sycl::ext::oneapi::level_zero::make<sycl::context>(sycl_devices, hContext,  sycl::ext::oneapi::level_zero::ownership::keep);
#endif
    for (int D=0; D< omp_get_num_devices(); D++)
      if (ompDeviceId2Context[D].ze_context == hContext)
        // This only work because the backend poiter is saved as a shared_pointer in SYCL context with Intel Implementation
        // https://github.com/intel/llvm/blob/ef33c57e48237c7d918f5dab7893554cecc001dd/sycl/source/backend/level_zero.cpp#L59
        // As far as I know this is not required by the SYCL2020 Spec
        ompDeviceId2Context[D].sycl_context = sycl_context;
  }
  return ompDeviceId2Context;
}

const struct syclDeviceInfo xomp_get_device_info(const int n) {
  return xomp_get_infos_devices()[n];
}

#include "openmc/event.h"
namespace openmc{

void sort_queue_SYCL(EventQueueItem* begin, EventQueueItem* end)
{
  static bool first {true};
  static sycl::queue sycl_command_queue;
  // Create a SYCL Queue for the device
  if (first) {
    int D = omp_get_default_device();
    sycl_command_queue = sycl::queue(xomp_get_device_info(D).sycl_context, xomp_get_device_info(D).sycl_device);
    first = false;
  }

  // Perform the sort on-device
  std::sort( oneapi::dpl::execution::make_device_policy(sycl_command_queue), begin, end);
}

} // end namespace openmc

#endif
