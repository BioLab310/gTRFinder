# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /data/run01/scw6eyx/miniconda3/envs/myenv/bin/cmake

# The command to remove a file.
RM = /data/run01/scw6eyx/miniconda3/envs/myenv/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /HOME/scw6eyx/run/cyf/gTRFinder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /HOME/scw6eyx/run/cyf/gTRFinder/build

# Include any dependencies generated for this target.
include CMakeFiles/gTRFinder.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/gTRFinder.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/gTRFinder.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gTRFinder.dir/flags.make

CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/DBGindex.cpp
CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o -MF CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/DBGindex.cpp

CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/DBGindex.cpp > CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.i

CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/DBGindex.cpp -o CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.s

CMakeFiles/gTRFinder.dir/src/basic.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/basic.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/basic.cpp
CMakeFiles/gTRFinder.dir/src/basic.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/gTRFinder.dir/src/basic.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/basic.cpp.o -MF CMakeFiles/gTRFinder.dir/src/basic.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/basic.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/basic.cpp

CMakeFiles/gTRFinder.dir/src/basic.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/basic.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/basic.cpp > CMakeFiles/gTRFinder.dir/src/basic.cpp.i

CMakeFiles/gTRFinder.dir/src/basic.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/basic.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/basic.cpp -o CMakeFiles/gTRFinder.dir/src/basic.cpp.s

CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/evaluate.cpp
CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o -MF CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/evaluate.cpp

CMakeFiles/gTRFinder.dir/src/evaluate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/evaluate.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/evaluate.cpp > CMakeFiles/gTRFinder.dir/src/evaluate.cpp.i

CMakeFiles/gTRFinder.dir/src/evaluate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/evaluate.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/evaluate.cpp -o CMakeFiles/gTRFinder.dir/src/evaluate.cpp.s

CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/generateSeq.cpp
CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o -MF CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/generateSeq.cpp

CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/generateSeq.cpp > CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.i

CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/generateSeq.cpp -o CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.s

CMakeFiles/gTRFinder.dir/src/hash.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/hash.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/hash.cpp
CMakeFiles/gTRFinder.dir/src/hash.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/gTRFinder.dir/src/hash.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/hash.cpp.o -MF CMakeFiles/gTRFinder.dir/src/hash.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/hash.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/hash.cpp

CMakeFiles/gTRFinder.dir/src/hash.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/hash.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/hash.cpp > CMakeFiles/gTRFinder.dir/src/hash.cpp.i

CMakeFiles/gTRFinder.dir/src/hash.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/hash.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/hash.cpp -o CMakeFiles/gTRFinder.dir/src/hash.cpp.s

CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/kmer_TRs_interval.cpp
CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o -MF CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/kmer_TRs_interval.cpp

CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/kmer_TRs_interval.cpp > CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.i

CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/kmer_TRs_interval.cpp -o CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.s

CMakeFiles/gTRFinder.dir/src/main.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/main.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/main.cpp
CMakeFiles/gTRFinder.dir/src/main.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/gTRFinder.dir/src/main.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/main.cpp.o -MF CMakeFiles/gTRFinder.dir/src/main.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/main.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/main.cpp

CMakeFiles/gTRFinder.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/main.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/main.cpp > CMakeFiles/gTRFinder.dir/src/main.cpp.i

CMakeFiles/gTRFinder.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/main.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/main.cpp -o CMakeFiles/gTRFinder.dir/src/main.cpp.s

CMakeFiles/gTRFinder.dir/src/match.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/match.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/match.cpp
CMakeFiles/gTRFinder.dir/src/match.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/gTRFinder.dir/src/match.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/match.cpp.o -MF CMakeFiles/gTRFinder.dir/src/match.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/match.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/match.cpp

CMakeFiles/gTRFinder.dir/src/match.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/match.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/match.cpp > CMakeFiles/gTRFinder.dir/src/match.cpp.i

CMakeFiles/gTRFinder.dir/src/match.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/match.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/match.cpp -o CMakeFiles/gTRFinder.dir/src/match.cpp.s

CMakeFiles/gTRFinder.dir/src/parallel.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/parallel.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/parallel.cpp
CMakeFiles/gTRFinder.dir/src/parallel.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/gTRFinder.dir/src/parallel.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/parallel.cpp.o -MF CMakeFiles/gTRFinder.dir/src/parallel.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/parallel.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/parallel.cpp

CMakeFiles/gTRFinder.dir/src/parallel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/parallel.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/parallel.cpp > CMakeFiles/gTRFinder.dir/src/parallel.cpp.i

CMakeFiles/gTRFinder.dir/src/parallel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/parallel.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/parallel.cpp -o CMakeFiles/gTRFinder.dir/src/parallel.cpp.s

CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/plot_data.cpp
CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o -MF CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/plot_data.cpp

CMakeFiles/gTRFinder.dir/src/plot_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/plot_data.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/plot_data.cpp > CMakeFiles/gTRFinder.dir/src/plot_data.cpp.i

CMakeFiles/gTRFinder.dir/src/plot_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/plot_data.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/plot_data.cpp -o CMakeFiles/gTRFinder.dir/src/plot_data.cpp.s

CMakeFiles/gTRFinder.dir/src/test.cpp.o: CMakeFiles/gTRFinder.dir/flags.make
CMakeFiles/gTRFinder.dir/src/test.cpp.o: /HOME/scw6eyx/run/cyf/gTRFinder/src/test.cpp
CMakeFiles/gTRFinder.dir/src/test.cpp.o: CMakeFiles/gTRFinder.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/gTRFinder.dir/src/test.cpp.o"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/gTRFinder.dir/src/test.cpp.o -MF CMakeFiles/gTRFinder.dir/src/test.cpp.o.d -o CMakeFiles/gTRFinder.dir/src/test.cpp.o -c /HOME/scw6eyx/run/cyf/gTRFinder/src/test.cpp

CMakeFiles/gTRFinder.dir/src/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/gTRFinder.dir/src/test.cpp.i"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /HOME/scw6eyx/run/cyf/gTRFinder/src/test.cpp > CMakeFiles/gTRFinder.dir/src/test.cpp.i

CMakeFiles/gTRFinder.dir/src/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/gTRFinder.dir/src/test.cpp.s"
	/HOME/scw6eyx/run/miniconda3/envs/myenv/bin/x86_64-conda-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /HOME/scw6eyx/run/cyf/gTRFinder/src/test.cpp -o CMakeFiles/gTRFinder.dir/src/test.cpp.s

# Object files for target gTRFinder
gTRFinder_OBJECTS = \
"CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/basic.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/hash.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/main.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/match.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/parallel.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o" \
"CMakeFiles/gTRFinder.dir/src/test.cpp.o"

# External object files for target gTRFinder
gTRFinder_EXTERNAL_OBJECTS =

gTRFinder: CMakeFiles/gTRFinder.dir/src/DBGindex.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/basic.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/evaluate.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/generateSeq.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/hash.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/kmer_TRs_interval.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/main.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/match.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/parallel.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/plot_data.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/src/test.cpp.o
gTRFinder: CMakeFiles/gTRFinder.dir/build.make
gTRFinder: /HOME/scw6eyx/run/miniconda3/envs/myenv/lib/libpython3.10.so
gTRFinder: CMakeFiles/gTRFinder.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX executable gTRFinder"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gTRFinder.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gTRFinder.dir/build: gTRFinder
.PHONY : CMakeFiles/gTRFinder.dir/build

CMakeFiles/gTRFinder.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gTRFinder.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gTRFinder.dir/clean

CMakeFiles/gTRFinder.dir/depend:
	cd /HOME/scw6eyx/run/cyf/gTRFinder/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /HOME/scw6eyx/run/cyf/gTRFinder /HOME/scw6eyx/run/cyf/gTRFinder /HOME/scw6eyx/run/cyf/gTRFinder/build /HOME/scw6eyx/run/cyf/gTRFinder/build /HOME/scw6eyx/run/cyf/gTRFinder/build/CMakeFiles/gTRFinder.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/gTRFinder.dir/depend

