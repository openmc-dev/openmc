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

# Utility rule file for ContinuousSubmit.

# Include any custom commands dependencies for this target.
include vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/compiler_depend.make

# Include the progress variables for this target.
include vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/progress.make

vendor/pugixml/CMakeFiles/ContinuousSubmit:
	cd /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml && /home/bamidele/.local/lib/python3.10/site-packages/cmake/data/bin/ctest -D ContinuousSubmit

ContinuousSubmit: vendor/pugixml/CMakeFiles/ContinuousSubmit
ContinuousSubmit: vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/build.make
.PHONY : ContinuousSubmit

# Rule to build all files generated by this target.
vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/build: ContinuousSubmit
.PHONY : vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/build

vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/clean:
	cd /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml && $(CMAKE_COMMAND) -P CMakeFiles/ContinuousSubmit.dir/cmake_clean.cmake
.PHONY : vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/clean

vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/depend:
	cd /mnt/c/Users/bamidele/openmc/build_release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/bamidele/openmc /mnt/c/Users/bamidele/openmc/vendor/pugixml /mnt/c/Users/bamidele/openmc/build_release /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : vendor/pugixml/CMakeFiles/ContinuousSubmit.dir/depend

