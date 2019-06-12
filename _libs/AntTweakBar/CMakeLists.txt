cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(AntTweakBar)

set(ATB_SOURCES 
	src/AntPerfTimer.h
	src/LoadOGL.cpp
	src/LoadOGL.h
	src/LoadOGLCore.cpp
	src/LoadOGLCore.h
	src/MiniGLFW.h
	src/MiniGLUT.h
	src/MiniSDL12.h
	src/MiniSDL13.h
	src/MiniSFML16.h
	src/TwBar.cpp
	src/TwBar.h
	src/TwColors.cpp
	src/TwColors.h
	src/TwEventGLFW.c
	src/TwEventGLUT.c
	src/TwEventSDL.c
	src/TwEventSDL12.c
	src/TwEventSDL13.c
	src/TwEventSFML.cpp
	src/TwFonts.cpp
	src/TwFonts.h
	src/TwGraph.h
	src/TwMgr.cpp
	src/TwMgr.h
	src/TwOpenGL.cpp
	src/TwOpenGL.h
	src/TwOpenGLCore.cpp
	src/TwOpenGLCore.h
	src/TwPrecomp.cpp
	src/TwPrecomp.h
)

set(ATB_HEADERS
    include/AntTweakBar.h
)
    
add_definitions(-DTW_EXPORTS -DTW_STATIC)

if (WIN32)
	REMOVE_DEFINITIONS(-DNOMINMAX)
    list(APPEND ATB_SOURCES src/TwEventWin.c)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd4005")
endif()

if (APPLE)
   include_directories ( /System/Library/Frameworks )
   find_library(APPKIT_LIBRARY AppKit)
   set(EXTRA_LIBS ${APPKIT_LIBRARY})
   set(CMAKE_CXX_FLAGS "-g -Wall -ObjC++ -fPIC -fno-strict-aliasing -D__PLACEMENT_NEW_INLINE")
   add_definitions(-D_MACOSX)   
elseif (UNIX)
    add_definitions(-D_UNIX -D__PLACEMENT_NEW_INLINE)
endif()

include_directories(include)

add_library(AntTweakBar STATIC ${ATB_HEADERS} ${ATB_SOURCES})
set_target_properties(AntTweakBar PROPERTIES FOLDER "_library")

if (APPLE)
    target_link_libraries(AntTweakBar) 
elseif (UNIX)
    target_link_libraries(AntTweakBar X11)
endif()