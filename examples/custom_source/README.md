# Building a Custom Source

To run this example, you first need to compile the custom source library, which
requires headers from OpenMC. A CMakeLists.txt file has been set up for you that
will search for OpenMC and build the custom library. To build the source
library, you can run:

    mkdir build && cd build
    OPENMC_ROOT=<path_to_openmc_install> cmake ..
    make

After this, you can build the model by running `python build_xml.py`. In the XML
files that are created, you should see a reference to build/libsource.so, the
custom source library that was built by CMake.
