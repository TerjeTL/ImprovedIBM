# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.24)

macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "-O2 /EHsc")

#	FETCHCONTENT METHOD
include(FetchContent)
FetchContent_Declare(
	catch
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v2.13.6
)
FetchContent_Declare(
	Eigen
	GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
	GIT_TAG           master
)


set(HDF5_BUILD_CPP_LIB ON)
set(BUILD_TESTING OFF)
set(HDF5_BUILD_EXAMPLES OFF)
FetchContent_Declare(
	HDF5
	URL		"https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.zip"
	OVERRIDE_FIND_PACKAGE
)

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
#FetchContent_MakeAvailable(catch)
FetchContent_MakeAvailable(Eigen)
# use target Eigen3::Eigen
#FetchContent_MakeAvailable(sdl2)
find_package(HDF5 COMPONENTS CXX REQUIRED)

FetchContent_GetProperties(highfive)
if(NOT highfive_POPULATED)
  FetchContent_Populate(highfive)
endif()

add_library(highfive INTERFACE)
target_include_directories(highfive INTERFACE ${highfive_SOURCE_DIR}/include)
#print_all_variables()

#	DOWNLOAD ALL THE SUBMODULES (git submodules)
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
# Update submodules as needed
	option(GIT_SUBMODULE "Check submodules during build" ON)
	if(GIT_SUBMODULE)
		message(STATUS "Submodule update")
		execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
						WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
						RESULT_VARIABLE GIT_SUBMOD_RESULT)
		if(NOT GIT_SUBMOD_RESULT EQUAL "0")
			message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT},
					please checkout submodules")
		endif()
	endif()
endif()

#	CHECK ALL THE SUBMODULES
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/SDL/CMakeLists.txt")
	message(FATAL_ERROR "The SDL submodule(s) was not downloaded! GIT_SUBMODULE was turned off or failed.")
endif()
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/SDL/CMakeLists.txt")
	message(FATAL_ERROR "The SDL2_ttf submodule(s) was not downloaded! GIT_SUBMODULE was turned off or failed.")
endif()
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/blaze/CMakeLists.txt")
	message(FATAL_ERROR "The Blaze submodule(s) was not downloaded! GIT_SUBMODULE was turned off or failed.")
endif()

project("ibm_application")

# Include sub-projects.
add_subdirectory("ibm_application")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/SDL")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/SDL_ttf")