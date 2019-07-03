cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(imgui)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fno-common -fPIC -Wall")
endif()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})
ADD_LIBRARY( imgui STATIC imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_widgets.cpp )
SET_PROPERTY(TARGET imgui PROPERTY FOLDER "_library")