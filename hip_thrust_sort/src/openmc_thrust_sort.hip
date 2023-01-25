#include<thrust/sort.h>

// This macro is used to enable the "__host__ __device__" attributes
// for the EventQueueItem comparator in event.h. If those attributes
// are missing, the code will compile and run, but it will use some sort
// of default comparator that does not have the desired effect.
#define COMPILE_CUDA_COMPARATOR

#include"openmc/event.h"

namespace openmc{

void device_sort_event_queue_item(EventQueueItem* begin, EventQueueItem* end)
{
  thrust::sort(thrust::device, begin, end);
}

}
