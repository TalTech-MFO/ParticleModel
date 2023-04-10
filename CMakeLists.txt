# Minimum version of CMake required
cmake_minimum_required(VERSION 3.12)

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

# Add the include directories
include_directories(${INC_DIR})

#=================================
# Biofouling implementation
#=================================
# Set the default implementation
# - simple
# - kooi
# - tsiaras
set(BIOFOULING_IMPL simple)
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

message(STATUS "Using netCDF4 include directory: ${NETCDF_INCLUDE_DIR}")
message(STATUS "Using netCDF4 libraries: ${NETCDF_LIBRARIES}")

# Append the netCDF4 include directories to the include target
include_directories(${NETCDF_INCLUDE_DIR})

#=================================
# Compiler flags
#=================================
# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()

# Add the debug flag in the debug configuration
if(CMAKE_BUILD_TYPE MATCHES Debug)
add_definitions(-DDEBUG)
endif()

# Set the compiler flags
set(CMAKE_Fortran_FLAGS "-cpp -Wall -ffree-line-length-none -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace -mcmodel=large")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSNAP_TO_BOUNDS -DSMAGORINSKY_INTERP_UV -DPARTICLE_SNAP_SEALVL -DIGNORE_BAD_PARTICLES")

# Set the optimization flags for the release configuration
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DUSE_OMP -fopenmp")

# Set the debugging flags for the debug configuration
set(CMAKE_Fortran_FLAGS_DEBUG "-Og -g")

# Add the executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the necessary libraries if needed
target_link_libraries(${PROJECT_NAME} ${LIBRARIES} ${NETCDF_LIBRARIES})
