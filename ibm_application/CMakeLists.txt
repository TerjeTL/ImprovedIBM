# CMakeList.txt : CMake project for gnuplot_application, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (ibm_application "src/ibm_application/ibm_application.cpp" "include/ibm_application/ibm_application.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ibm_application PROPERTY CXX_STANDARD 20)
endif()

target_include_directories(ibm_application
	PUBLIC "${PROJECT_SOURCE_DIR}/ibm_application/include/"
)

set(SDL_LIBRARIES ${SDL_LIBRARIES} SDL2main SDL2-static)
target_link_libraries(ibm_application
	PRIVATE ${SDL_LIBRARIES}
)
# TODO: Add tests and install targets if needed.
