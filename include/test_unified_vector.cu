#include "openmc/settings.h"
#include "openmc/vector.h"
#include <iostream>
#include <type_traits>

using namespace openmc;
bool gpu_mode = true;

__global__ void set_values(vector<double> data)
{
  for (unsigned i = 0; i < data.size(); i++)
    data[i] = i + 1;
}

int main()
{
  vector<double> x(5, 1);
  for (unsigned i = 0; i < x.size(); i++)
    std::cout << x[i] << std::endl;
  std::cout << "set through cuda" << std::endl;
  set_values<<<1, 1>>>(x);
  cudaDeviceSynchronize();
  for (unsigned i = 0; i < x.size(); i++)
    std::cout << x[i] << std::endl;

  std::cout << "push back on host, then reset" << std::endl;
  x.push_back(123);
  set_values<<<1, 1>>>(x);
  cudaDeviceSynchronize();
  for (unsigned i = 0; i < x.size(); i++)
    std::cout << x[i] << std::endl;
}
