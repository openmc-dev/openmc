# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/bamidele/.local/lib/python3.10/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/bamidele/.local/lib/python3.10/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/bamidele/openmc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/bamidele/openmc/build_release

# Include any dependencies generated for this target.
include tests/cpp_unit_tests/CMakeFiles/test_tally.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/cpp_unit_tests/CMakeFiles/test_tally.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/cpp_unit_tests/CMakeFiles/test_tally.dir/progress.make

# Include the compile flags for this target's objects.
include tests/cpp_unit_tests/CMakeFiles/test_tally.dir/flags.make

tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o: tests/cpp_unit_tests/CMakeFiles/test_tally.dir/flags.make
tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o: /mnt/c/Users/bamidele/openmc/tests/cpp_unit_tests/test_tally.cpp
tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o: tests/cpp_unit_tests/CMakeFiles/test_tally.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/mnt/c/Users/bamidele/openmc/build_release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o"
	cd /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o -MF CMakeFiles/test_tally.dir/test_tally.cpp.o.d -o CMakeFiles/test_tally.dir/test_tally.cpp.o -c /mnt/c/Users/bamidele/openmc/tests/cpp_unit_tests/test_tally.cpp

tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/test_tally.dir/test_tally.cpp.i"
	cd /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/bamidele/openmc/tests/cpp_unit_tests/test_tally.cpp > CMakeFiles/test_tally.dir/test_tally.cpp.i

tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/test_tally.dir/test_tally.cpp.s"
	cd /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/bamidele/openmc/tests/cpp_unit_tests/test_tally.cpp -o CMakeFiles/test_tally.dir/test_tally.cpp.s

# Object files for target test_tally
test_tally_OBJECTS = \
"CMakeFiles/test_tally.dir/test_tally.cpp.o"

# External object files for target test_tally
test_tally_EXTERNAL_OBJECTS =

bin/test_tally: tests/cpp_unit_tests/CMakeFiles/test_tally.dir/test_tally.cpp.o
bin/test_tally: tests/cpp_unit_tests/CMakeFiles/test_tally.dir/build.make
bin/test_tally: lib/libCatch2Main.a
bin/test_tally: lib/libopenmc.so
bin/test_tally: lib/libCatch2.a
bin/test_tally: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
bin/test_tally: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.so
bin/test_tally: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
bin/test_tally: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.so
bin/test_tally: lib/libfmt.a
bin/test_tally: lib/libpugixml.a
bin/test_tally: /usr/lib/x86_64-linux-gnu/libpng.so
bin/test_tally: /usr/lib/x86_64-linux-gnu/libz.so
bin/test_tally: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
bin/test_tally: /usr/lib/x86_64-linux-gnu/libpthread.so
bin/test_tally: tests/cpp_unit_tests/CMakeFiles/test_tally.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/mnt/c/Users/bamidele/openmc/build_release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/test_tally"
	cd /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_tally.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/cpp_unit_tests/CMakeFiles/test_tally.dir/build: bin/test_tally
.PHONY : tests/cpp_unit_tests/CMakeFiles/test_tally.dir/build

tests/cpp_unit_tests/CMakeFiles/test_tally.dir/clean:
	cd /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests && $(CMAKE_COMMAND) -P CMakeFiles/test_tally.dir/cmake_clean.cmake
.PHONY : tests/cpp_unit_tests/CMakeFiles/test_tally.dir/clean

tests/cpp_unit_tests/CMakeFiles/test_tally.dir/depend:
	cd /mnt/c/Users/bamidele/openmc/build_release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/bamidele/openmc /mnt/c/Users/bamidele/openmc/tests/cpp_unit_tests /mnt/c/Users/bamidele/openmc/build_release /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests /mnt/c/Users/bamidele/openmc/build_release/tests/cpp_unit_tests/CMakeFiles/test_tally.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tests/cpp_unit_tests/CMakeFiles/test_tally.dir/depend

