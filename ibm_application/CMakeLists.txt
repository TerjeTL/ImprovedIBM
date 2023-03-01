# CMakeList.txt : CMake project for gnuplot_application, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.24)

find_package(OpenMP REQUIRED)

# Add source to this project's executable.
add_executable(ibm_application "src/ibm_application/ibm_application.cpp")

target_sources(ibm_application PRIVATE
	"src/ibm_application/DataViewer.cpp"
	"src/ibm_application/GeometrySDF.cpp"
	"src/ibm_application/CartGrid.cpp"
	"src/ibm_application/Schemes.cpp"
	"src/ibm_application/Solver.cpp"
	"src/ibm_application/DataExporter.cpp"
	"src/data_viewer/Shader.cpp"
	"src/data_viewer/Camera.cpp"
	"src/data_viewer/Mesh.cpp"
	"src/data_viewer/SolutionModel.cpp"
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ibm_application PROPERTY CXX_STANDARD 17)
endif()

target_compile_definitions(ibm_application PRIVATE _USE_MATH_DEFINES)

target_include_directories(ibm_application
	PUBLIC "${PROJECT_SOURCE_DIR}/ibm_application/include/"
#	PUBLIC "${CMAKE_BINARY_DIR}/_deps/hdf5-src/src/H5FDsubfiling"
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/blaze/cmake")
#include(Blaze_Import)
#Blaze_Import(
#	QUIET
#	BLAS off
#	LAPACK off
#	THREADING off
#	CACHE_SIZE auto
#	VECTORIZATION off
#	STORAGE_ORDER rowMajor
#)

set(SDL_LIBRARIES ${SDL_LIBRARIES} SDL2main SDL2 SDL2_ttf)

target_link_libraries(ibm_application
	PUBLIC OpenMP::OpenMP_CXX
	PRIVATE ${SDL_LIBRARIES}
#	PRIVATE Blaze
	PRIVATE Eigen3::Eigen
	PRIVATE hdf5_data_export
	imgui
)

#print_all_variables()
# TODO: Add tests and install targets if needed.
