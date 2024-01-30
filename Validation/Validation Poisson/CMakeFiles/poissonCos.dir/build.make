# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Validation poisson"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Validation poisson"

# Include any dependencies generated for this target.
include CMakeFiles/poissonCos.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poissonCos.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poissonCos.dir/flags.make

CMakeFiles/poissonCos.dir/poissonCos.cpp.o: CMakeFiles/poissonCos.dir/flags.make
CMakeFiles/poissonCos.dir/poissonCos.cpp.o: poissonCos.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Validation poisson/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/poissonCos.dir/poissonCos.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonCos.dir/poissonCos.cpp.o -c "/Validation poisson/poissonCos.cpp"

CMakeFiles/poissonCos.dir/poissonCos.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonCos.dir/poissonCos.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Validation poisson/poissonCos.cpp" > CMakeFiles/poissonCos.dir/poissonCos.cpp.i

CMakeFiles/poissonCos.dir/poissonCos.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonCos.dir/poissonCos.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Validation poisson/poissonCos.cpp" -o CMakeFiles/poissonCos.dir/poissonCos.cpp.s

CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o: CMakeFiles/poissonCos.dir/flags.make
CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o: prediction-projection/Prediction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Validation poisson/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o -c "/Validation poisson/prediction-projection/Prediction.cpp"

CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Validation poisson/prediction-projection/Prediction.cpp" > CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.i

CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Validation poisson/prediction-projection/Prediction.cpp" -o CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.s

CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o: CMakeFiles/poissonCos.dir/flags.make
CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o: prediction-projection/Projection.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Validation poisson/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o -c "/Validation poisson/prediction-projection/Projection.cpp"

CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Validation poisson/prediction-projection/Projection.cpp" > CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.i

CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Validation poisson/prediction-projection/Projection.cpp" -o CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.s

# Object files for target poissonCos
poissonCos_OBJECTS = \
"CMakeFiles/poissonCos.dir/poissonCos.cpp.o" \
"CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o" \
"CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o"

# External object files for target poissonCos
poissonCos_EXTERNAL_OBJECTS =

poissonCos: CMakeFiles/poissonCos.dir/poissonCos.cpp.o
poissonCos: CMakeFiles/poissonCos.dir/prediction-projection/Prediction.cpp.o
poissonCos: CMakeFiles/poissonCos.dir/prediction-projection/Projection.cpp.o
poissonCos: CMakeFiles/poissonCos.dir/build.make
poissonCos: /usr/local/lib/libneos.so
poissonCos: /usr/local/lib/libbitpit_MPI.so
poissonCos: /usr/local/petsc/lib/libpetsc.so
poissonCos: /usr/lib/x86_64-linux-gnu/libxml2.so
poissonCos: /usr/lib/x86_64-linux-gnu/liblapacke.so
poissonCos: /usr/lib/x86_64-linux-gnu/liblapack.so
poissonCos: /usr/lib/x86_64-linux-gnu/libblas.so
poissonCos: /usr/local/petsc/lib/libpetsc.so
poissonCos: /usr/lib/x86_64-linux-gnu/libxml2.so
poissonCos: /usr/lib/x86_64-linux-gnu/liblapacke.so
poissonCos: /usr/lib/x86_64-linux-gnu/liblapack.so
poissonCos: /usr/lib/x86_64-linux-gnu/libblas.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libstdc++.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libgcc_s.so
poissonCos: /usr/lib/x86_64-linux-gnu/libpthread.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libstdc++.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so
poissonCos: /usr/lib/gcc/x86_64-linux-gnu/9/libgcc_s.so
poissonCos: /usr/lib/x86_64-linux-gnu/libpthread.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
poissonCos: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
poissonCos: CMakeFiles/poissonCos.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Validation poisson/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable poissonCos"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poissonCos.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poissonCos.dir/build: poissonCos

.PHONY : CMakeFiles/poissonCos.dir/build

CMakeFiles/poissonCos.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poissonCos.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poissonCos.dir/clean

CMakeFiles/poissonCos.dir/depend:
	cd "/Validation poisson" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Validation poisson" "/Validation poisson" "/Validation poisson" "/Validation poisson" "/Validation poisson/CMakeFiles/poissonCos.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/poissonCos.dir/depend
