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

# Utility rule file for NightlySubmit.

# Include any custom commands dependencies for this target.
include vendor/pugixml/CMakeFiles/NightlySubmit.dir/compiler_depend.make

# Include the progress variables for this target.
include vendor/pugixml/CMakeFiles/NightlySubmit.dir/progress.make

vendor/pugixml/CMakeFiles/NightlySubmit:
	cd /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml && /home/bamidele/.local/lib/python3.10/site-packages/cmake/data/bin/ctest -D NightlySubmit

NightlySubmit: vendor/pugixml/CMakeFiles/NightlySubmit
NightlySubmit: vendor/pugixml/CMakeFiles/NightlySubmit.dir/build.make
.PHONY : NightlySubmit

# Rule to build all files generated by this target.
vendor/pugixml/CMakeFiles/NightlySubmit.dir/build: NightlySubmit
.PHONY : vendor/pugixml/CMakeFiles/NightlySubmit.dir/build

vendor/pugixml/CMakeFiles/NightlySubmit.dir/clean:
	cd /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml && $(CMAKE_COMMAND) -P CMakeFiles/NightlySubmit.dir/cmake_clean.cmake
.PHONY : vendor/pugixml/CMakeFiles/NightlySubmit.dir/clean

vendor/pugixml/CMakeFiles/NightlySubmit.dir/depend:
	cd /mnt/c/Users/bamidele/openmc/build_release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/bamidele/openmc /mnt/c/Users/bamidele/openmc/vendor/pugixml /mnt/c/Users/bamidele/openmc/build_release /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml /mnt/c/Users/bamidele/openmc/build_release/vendor/pugixml/CMakeFiles/NightlySubmit.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : vendor/pugixml/CMakeFiles/NightlySubmit.dir/depend

