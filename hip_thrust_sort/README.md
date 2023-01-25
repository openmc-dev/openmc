This subdirectory contains code and build routines for
performing on-device sorting via the CUDA Thrust library.
The use of a separate subdirectory is necessary as OpenMP
compilers do not know how to link link to Thrust. To avoid
polluting the main CMakeLists.txt for OpenMC, we have used
a subdirectory instead to isolate the build for just the 
thrust portion and then to link to the library that is 
generated.
