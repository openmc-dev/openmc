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
custom source library that was built by CMake. The model is also set up with a
mesh tally of the flux, so once you run `openmc`, you will get a statepoint file
with the tally results in it. Running `python show_flux.py` will pull in the
results from the statepoint file and display them. If all worked well, you
should see a ring "imprint" as well as a higher flux to the right side (since
the custom source has all particles moving in the positive x direction).
