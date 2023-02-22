#ifdef SYCL_SORT
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <omp.h>

#include "openmc/event.h"

namespace openmc{

void sort_queue_SYCL(EventQueueItem* begin, EventQueueItem* end)
{
  std::sort( oneapi::dpl::execution::dpcpp_default, begin, end);
}

} // end namespace openmc

#endif
