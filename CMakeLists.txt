# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.14)

#	FETCHCONTENT METHOD
include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.6
)
#FetchContent_Declare(
#  sdl2
#  GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
#  GIT_TAG        release-2.24.x
#)

# CMake 3.14+
FetchContent_MakeAvailable(catch)
#FetchContent_MakeAvailable(sdl2)

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
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/blaze/CMakeLists.txt")
	message(FATAL_ERROR "The Blaze submodule(s) was not downloaded! GIT_SUBMODULE was turned off or failed.")
endif()

project("ibm_application")

# Include sub-projects.
add_subdirectory("ibm_application")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/SDL")