cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(tinyxml)

SET(SRC_TINYXML
	tinystr.cpp
	tinyxml.cpp
	tinyxmlerror.cpp
	tinyxmlparser.cpp
	tinystr.h
	tinyxml.h
)

add_definitions(-D_LIB)
add_definitions(-DTIXML_USE_STL)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fno-common -fPIC -Wall")
endif()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
ADD_LIBRARY( tinyxml STATIC ${SRC_TINYXML} )
set_target_properties(tinyxml PROPERTIES FOLDER "_library")