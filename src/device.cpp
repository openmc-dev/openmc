#include "openmc/device.h"

#include "openmc/nuclide.h"
#include "openmc/output.h"

namespace openmc {

void map_to_device()
{
  write_message("Copying data to device...", 7);
  for (auto& nuc : data::nuclides) {
    for (auto& rx : nuc->reactions_) {
      for (auto& product : rx->products_) {
        for (auto& d : product.distribution_) {
          #pragma omp target enter data map(to: d)
          d.copy_to_device();
        }
      }
    }
  }
}

void release_from_device()
{
  write_message("Releasing data from device...", 7);
  for (auto& nuc : data::nuclides) {
    for (auto& rx : nuc->reactions_) {
      for (auto& product : rx->products_) {
        for (auto& d : product.distribution_) {
          d.release_device();
          #pragma omp target exit data map(release: d)
        }
      }
    }
  }
}

} // end namespace openmc
