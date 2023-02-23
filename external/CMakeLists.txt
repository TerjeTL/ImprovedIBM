cmake_minimum_required (VERSION 3.24)

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

add_subdirectory("${CMAKE_SOURCE_DIR}/external/SDL")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/SDL_ttf")

set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
set(SDL_DIR ${CMAKE_SOURCE_DIR}/external/SDL)
set(GLM_DIR ${CMAKE_SOURCE_DIR}/external/glm)
set(GLAD_DIR ${CMAKE_SOURCE_DIR}/external/glad)

add_library(imgui STATIC)
target_sources( imgui
                PRIVATE
                    ${IMGUI_DIR}/imgui_demo.cpp
                    ${IMGUI_DIR}/imgui_draw.cpp
                    ${IMGUI_DIR}/imgui_tables.cpp
                    ${IMGUI_DIR}/imgui_widgets.cpp
                    ${IMGUI_DIR}/imgui.cpp

                PRIVATE
                    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
                    ${IMGUI_DIR}/backends/imgui_impl_sdl2.cpp
                )

target_include_directories( imgui
							PUBLIC ${IMGUI_DIR}
							PUBLIC ${IMGUI_DIR}/backends
							PUBLIC ${SDL_DIR}/include
							PUBLIC ${GLM_DIR}/glm
)

find_package(OpenGL REQUIRED)
target_link_libraries(imgui PUBLIC ${OPENGL_LIBRARIES})

add_library(glad STATIC ${GLAD_DIR}/include/glad/glad.h ${GLAD_DIR}/src/glad.c)
target_include_directories(glad PUBLIC ${GLAD_DIR}/include/)

target_link_libraries(imgui PUBLIC SDL2 ${CMAKE_DL_LIBS} glad)