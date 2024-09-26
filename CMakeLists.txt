# Minimum version of CMake required
cmake_minimum_required(VERSION 3.12)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the name of the project
project(ParticleModel Fortran)

#=================================
# Source paths
#=================================
# Set the path to the source code directory
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/code)

# Set the path to the include directory
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add the source files
file(GLOB SOURCES ${SRC_DIR}/*.f90)
option(USE_VECTORIZED_POSTPROCESSING "Use the vectorized postprocessing implementation" OFF)
if (USE_VECTORIZED_POSTPROCESSING)
    message(STATUS "Using vectorized postprocessing implementation")
    set(PP_DIR ${SRC_DIR}/postprocessing_vectorized)
    # Add a preprocessor definition
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUSE_VECTORIZED_POSTPROCESSING")
else()
    message(STATUS "Using non-vectorized postprocessing implementation")
    set(PP_DIR ${SRC_DIR}/postprocessing)
endif()
file(GLOB PP_SOURCES ${PP_DIR}/*.f90)
list(APPEND SOURCES ${PP_SOURCES})

# Add the include directories
include_directories(${INC_DIR})

#=================================
# Biofouling implementation
#=================================
# Set the default implementation
# - chl_proxy
set(BIOFOULING_IMPL chl_proxy)
message(STATUS "Using biofouling implementation: ${BIOFOULING_IMPL}")
file(GLOB IMPL_SOURCES ${SRC_DIR}/biofouling/${BIOFOULING_IMPL}/*.f90)
message(STATUS "Using biofouling implementation source files: ${IMPL_SOURCES}")
list(APPEND SOURCES ${IMPL_SOURCES})

#=================================
# External libraries
#=================================
# Find the netCDF4 include directories and libraries
execute_process(COMMAND nf-config --includedir OUTPUT_VARIABLE NETCDF_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --flibs OUTPUT_VARIABLE NETCDF_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "netCDF4 include directory: ${NETCDF_INCLUDE_DIR}")
message(STATUS "netCDF4 libraries: ${NETCDF_LIBRARIES}")

# Append the netCDF4 include directories to the include target
include_directories(${NETCDF_INCLUDE_DIR})

#=================================
# Compiler flags
#=================================
# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set the compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -Wall -ffree-line-length-none -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large -march=native")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSNAP_TO_BOUNDS -DPARTICLE_SNAP_SEALVL -DIGNORE_BAD_PARTICLES -DBIOFOULING_SIMPLE -DNO_DIFFUSE_VERTICAL")

# Set the optimization flags for the release configuration
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -DUSE_OMP -fopenmp")

# Set the debugging flags for the debug configuration
set(CMAKE_Fortran_FLAGS_DEBUG "-O3 -g -DDEBUG")

# Set the profile flag
option(BUILD_FOR_PROFILING "Build with -pg and -no-pie flags" OFF)
if (BUILD_FOR_PROFILING)
    message(STATUS "Configuring with -pg and -no-pie")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg -no-pie")
endif()

# Add the executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the necessary libraries if needed
target_link_libraries(${PROJECT_NAME} ${LIBRARIES} ${NETCDF_LIBRARIES})
