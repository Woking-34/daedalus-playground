cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(stb_image)

SET(SRC_STB_IMAGE
	stb_image.c
	stb_image.h
)

add_definitions(-DSTBI_NO_STDIO=1)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fno-common -fPIC -Wall")
endif()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
ADD_LIBRARY( stbi STATIC ${SRC_STB_IMAGE} )
SET_PROPERTY(TARGET stbi PROPERTY FOLDER "_library")