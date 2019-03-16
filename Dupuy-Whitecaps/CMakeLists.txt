SET(SRC_DUPUY_WHITECAPS
	main.cpp
	
	vec4.h
	mat4.h
	
	Program.h
	Program.cpp
)

add_executable(dupuy_whitecaps ${SRC_DUPUY_WHITECAPS})

IF(WIN32)
target_link_libraries(dupuy_whitecaps freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(dupuy_whitecaps glut GL)
ENDIF(UNIX AND NOT APPLE)

target_link_libraries(dupuy_whitecaps glew)
target_link_libraries(dupuy_whitecaps AntTweakBar)

target_compile_definitions(dupuy_whitecaps PRIVATE SAMPLE_NAME="Dupuy-Whitecaps")
