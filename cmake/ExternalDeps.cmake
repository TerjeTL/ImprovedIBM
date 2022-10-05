cmake_minimum_required (VERSION 3.24)

#	FETCHCONTENT METHOD
include(FetchContent)
# Fetch catch test
FetchContent_Declare(
	catch
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v2.13.6
)
# Fetch Eigen linear algebra library
FetchContent_Declare(
	Eigen
	GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
	GIT_TAG           master
)

# Fetch HDF5 library
set(HDF5_BUILD_CPP_LIB ON)
set(BUILD_TESTING OFF)
set(HDF5_BUILD_EXAMPLES OFF)
FetchContent_Declare(
	HDF5
	URL		"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.zip"
	OVERRIDE_FIND_PACKAGE
)

# Fetch HighFive header only HDF5 cpp high-level library
Fetchcontent_Declare(
	highfive
	URL "https://github.com/BlueBrain/HighFive/archive/refs/tags/v2.4.1.zip"
)

#FetchContent_Declare(
#  sdl2
#  GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
#  GIT_TAG        release-2.24.x
#)

# CMake 3.14+ 
FetchContent_MakeAvailable(Eigen)


#	DATA EXPORT LIBRARIES
find_package(HDF5 COMPONENTS CXX REQUIRED)
FetchContent_GetProperties(highfive)
if(NOT highfive_POPULATED)
  FetchContent_Populate(highfive)
endif()

add_library(highfive INTERFACE)
target_include_directories(highfive INTERFACE ${highfive_SOURCE_DIR}/include)

# Interface target for data export libraries
add_library(hdf5_data_export INTERFACE)
target_include_directories(hdf5_data_export
	INTERFACE "${CMAKE_BINARY_DIR}/_deps/hdf5-src/src/H5FDsubfiling"
)
target_link_libraries(hdf5_data_export INTERFACE
	hdf5-shared
	hdf5_cpp-shared
	highfive
)