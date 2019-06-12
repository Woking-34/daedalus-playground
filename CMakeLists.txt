cmake_minimum_required(VERSION 2.6)
project(Daedalus-Playground)

set( CMAKE_VERBOSE_MAKEFILE ON )
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "_cmake")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/daedalus-bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/daedalus-lib)

IF(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
ENDIF()

string(TOUPPER ${CMAKE_CXX_COMPILER_ID} CXX_COMPILER_ID_UPPER)
IF(${CXX_COMPILER_ID_UPPER} STREQUAL MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP")
ENDIF()

add_definitions(-DTIXML_USE_STL)
add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

IF(WIN32)
	add_definitions(-DWIN32_LEAN_AND_MEAN)
	add_definitions(-DNOMINMAX)
	add_definitions(-D_USE_MATH_DEFINES)
	
	add_definitions(-DGLEW_STATIC)
	add_definitions(-DFREEGLUT_STATIC)
	add_subdirectory(_libs/freeglut "${CMAKE_CURRENT_BINARY_DIR}/_libs/freeglut")
	
	include_directories(_libs/freeglut)
ENDIF(WIN32)

add_subdirectory(_libs/glew "${CMAKE_CURRENT_BINARY_DIR}/_libs/glew")
add_subdirectory(_libs/AntTweakBar "${CMAKE_CURRENT_BINARY_DIR}/_libs/AntTweakBar")
add_subdirectory(_libs/tinyxml "${CMAKE_CURRENT_BINARY_DIR}/_libs/tinyxml")

include_directories(
	.

	_libs/glew
	_libs/AntTweakBar/include
	_libs/tinyxml
	_libs/glm

	_libs/eigen
	_libs/eigen/Eigen
)

find_package(CUDA)
IF(CUDA_FOUND)
	MESSAGE(STATUS "CUDA_TOOLKIT_ROOT_DIR - ${CUDA_TOOLKIT_ROOT_DIR}")

	SET(CUDA_VERBOSE_BUILD ON)
	LIST(APPEND CUDA_NVCC_FLAGS "--use_fast_math")
	LIST(APPEND CUDA_NVCC_FLAGS "--ptxas-options=-v")

	include_directories(${CUDA_INCLUDE_DIRS})
ELSE()
	MESSAGE(STATUS "CUDA NOT FOUND!")
ENDIF()

add_subdirectory(Bruneton-Ocean)
add_subdirectory(Dupuy-Whitecaps)
add_subdirectory(FastMassSpring)

IF(CUDA_FOUND AND WIN32)
	add_subdirectory(HuaminWang-Chebyshev)
	add_subdirectory(HuaminWang-DescentMethods)
ENDIF()
