#include<thrust/sort.h>

#define COMPILE_CUDA_COMPARATOR

#include"openmc/event.h"

namespace openmc{

  /*
struct EventQueueItem{
  int idx;         //!< particle index in event-based particle buffer
  //int material;    //!< material that particle is in
  float E;            //!< particle energy

  // Compare by particle type, then by material type (4.5% fuel/7.0% fuel/cladding/etc),
  // then by energy.
  // TODO: Currently in OpenMC, the material ID corresponds not only to a general
  // type, but also specific isotopic densities. Ideally we would
  // like to be able to just sort by general material type, regardless of densities.
  // A more general material type ID may be added in the future, in which case we
  // can update the material field of this struct to contain the more general id.
  __host__ __device__
  bool operator<(const EventQueueItem& rhs) const
  {
    return E < rhs.E;
  }

  // This is needed by the implementation of parallel quicksort
  __host__ __device__
  bool operator>(const EventQueueItem& rhs) const
  {
    return E > rhs.E;
  }
};
*/

void device_sort_event_queue_item(EventQueueItem* begin, EventQueueItem* end)
{
  thrust::sort(thrust::device, begin, end);
}

}
