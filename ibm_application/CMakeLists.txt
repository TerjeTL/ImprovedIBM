# CMakeList.txt : CMake project for gnuplot_application, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (ibm_application "src/ibm_application/ibm_application.cpp" "include/ibm_application/ibm_application.h" "src/ibm_application/SDLGraphics.cpp" "include/ibm_application/SDLGraphics.h" "src/ibm_application/GeometrySDF.cpp" "include/ibm_application/GeometrySDF.h" "include/ibm_application/CartGrid.h" "src/ibm_application/CartGrid.cpp" "include/ibm_application/Schemes.h" "src/ibm_application/Schemes.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ibm_application PROPERTY CXX_STANDARD 17)
endif()


target_include_directories(ibm_application
	PUBLIC "${PROJECT_SOURCE_DIR}/ibm_application/include/"
	PUBLIC "${CMAKE_SOURCE_DIR}/external/SDL/include/"
	PUBLIC "${CMAKE_SOURCE_DIR}/external/SDL_ttf/"
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/blaze/cmake")
include(Blaze_Import)
Blaze_Import(
	QUIET
	BLAS off
	LAPACK off
	THREADING off
	CACHE_SIZE auto
	VECTORIZATION off
	STORAGE_ORDER rowMajor
)

set(SDL_LIBRARIES ${SDL_LIBRARIES} SDL2main SDL2 SDL2_ttf)
target_link_libraries(ibm_application
	PRIVATE ${SDL_LIBRARIES}
	PRIVATE Blaze
	PRIVATE Eigen3::Eigen
)
# TODO: Add tests and install targets if needed.
