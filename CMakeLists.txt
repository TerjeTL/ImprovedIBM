# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.24)

#	CMAKE Defines
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "-O2 /EHsc /openmp")

#	CMAKE Utilities
include("${CMAKE_SOURCE_DIR}/cmake/macros.cmake")

#	CMAKE Files
include("${CMAKE_SOURCE_DIR}/cmake/ExternalDeps.cmake")
include("${CMAKE_SOURCE_DIR}/external/CMakeLists.txt")

#	CMAKE Project
project("ibm_application")

#	SUBDIRECTORIES
add_subdirectory("ibm_application")