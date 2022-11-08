////////////////////////////////////////////////////////////////////////////////////
// Parallel Quicksort Sorting Algorithm Using Task-Based OpenMP Threading 
////////////////////////////////////////////////////////////////////////////////////
//
// These algorithms are based on the parallel quicksort implementation by
// Eduard Lopez published at https://github.com/eduardlopez/quicksort-parallel
//
// Eduard's original version was for an integer type quicksort, but I have modified
// it to be templated. Additionally, I have modified the
// optimal chunk sizes and restricted the number of threads for the array sizing
// that XSBench will be using by default.
//
// Eduard's original implementation carries the following license, which applies to
// the following functions only:
//
//	void quickSort_parallel_internal()
//  void quickSort_parallel()
//
// The MIT License (MIT)
//
// Copyright (c) 2016 Eduard López
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////////

#include<algorithm>

template <class T>
void quickSort_parallel_internal(T* arr, int left, int right, int cutoff)
{
  int i = left;
  int j = right;
  T pivot = arr[(left + right) / 2];

  while (i <= j) {
    while (arr[i] < pivot)
      i++;
    while (arr[j] > pivot)
      j--;
    if (i <= j) {
      std::swap(arr[i], arr[j]);
      i++;
      j--;
    }
  }

  if ((right-left) < cutoff) {
    if (left < j)  { quickSort_parallel_internal(arr, left, j,  cutoff); }
    if (i < right) { quickSort_parallel_internal(arr, i, right, cutoff); }
  } else {
    #pragma omp task
    { quickSort_parallel_internal(arr, left, j, cutoff); }
    #pragma omp task
    { quickSort_parallel_internal(arr, i, right, cutoff); }
  }

}

template <class T>
void quickSort_parallel(T* arr, int lenArray)
{
  // Set minumum problem size to still spawn threads for
  int cutoff = 1000;

  // For OpenMC's use case, we do not see performance gains past 32 threads.
  // So, we will limit the number of threads to 32. 
  #ifdef _OPENMP
  int	numThreads = std::min({32,omp_get_num_procs()});
  #else
  int	numThreads = 1;
  #endif

  #pragma omp parallel num_threads(numThreads)
  {
    #pragma omp single nowait
    { quickSort_parallel_internal(arr, 0, lenArray-1, cutoff); }
  }
}
