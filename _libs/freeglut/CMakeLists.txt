cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(freeglut)

SET(SRC_FREEGLUT
	freeglut_callbacks.c
	freeglut_cursor.c
	freeglut_display.c
	freeglut_ext.c
	freeglut_font.c
	freeglut_font_data.c
	freeglut_gamemode.c
	freeglut_geometry.c
	freeglut_glutfont_definitions.c
	freeglut_init.c
	freeglut_input_devices.c
	freeglut_joystick.c
	freeglut_main.c
	freeglut_menu.c
	freeglut_misc.c
	freeglut_overlay.c
	freeglut_spaceball.c
	freeglut_state.c
	freeglut_stroke_mono_roman.c
	freeglut_stroke_roman.c
	freeglut_structure.c
	freeglut_teapot.c
	freeglut_videoresize.c
	freeglut_window.c
	freeglut_xinput.c
)

if(MSVC)
  set(DEF_FILE freeglutdll.def)
endif()

ADD_DEFINITIONS(-DFREEGLUT_EXPORTS)
ADD_DEFINITIONS(-DFREEGLUT_STATIC)
ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)

if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fno-common -fPIC -Wall")
endif()

include_directories(${CMAKE_CURRENT_LIST_DIR})
add_library( freeglut STATIC ${SRC_FREEGLUT} ${DEF_FILE} )
set_target_properties(freeglut PROPERTIES FOLDER "External Dependencies")