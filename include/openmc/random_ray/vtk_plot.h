#ifndef OPENMC_RANDOM_RAY_VTK_PLOT_H
#define OPENMC_RANDOM_RAY_VTK_PLOT_H

#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {
void plot_3D_vtk();

// Returns the inputted value with its
// endianness reversed. This is useful
// for conforming to the paraview VTK
// binary file format.
template<typename T>
T flip_endianness(T in)
{
  char* orig = reinterpret_cast<char*>(&in);
  char swapper[sizeof(T)];
  for (int i = 0; i < sizeof(T); i++) {
    swapper[i] = orig[sizeof(T) - i - 1];
  }
  T out = *reinterpret_cast<T*>(&swapper);
  return out;
}

} // namespace openmc

#endif // OPEMMC_RANMDOM_RAY_VTK_PLOT_H
