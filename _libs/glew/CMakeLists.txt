cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(glew)

SET(SRC_GLEW
	glew.c
)

add_definitions(-DGLEW_BUILD)
add_definitions(-DGLEW_NO_GLU)
add_definitions(-DGLEW_STATIC)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fno-common -fPIC -Wall")
endif()

include_directories(${CMAKE_CURRENT_LIST_DIR})
add_library( glew STATIC ${SRC_GLEW} )
set_target_properties(glew PROPERTIES FOLDER "External Dependencies")