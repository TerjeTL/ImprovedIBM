# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.14)

include(FetchContent)

FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.6
)
FetchContent_Declare(
  sdl2
  GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
  GIT_TAG        release-2.24.x
)

# CMake 3.14+
FetchContent_MakeAvailable(catch)
FetchContent_MakeAvailable(sdl2)

project ("ibm_application")

# Include sub-projects.
add_subdirectory ("ibm_application")

