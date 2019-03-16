SET(SRC_BRUNETON_OCEAN
	main.cpp
	
	vec4.h
	mat4.h
	
	Program.h
	Program.cpp
)

add_executable(bruneton_ocean ${SRC_BRUNETON_OCEAN})

IF(WIN32)
target_link_libraries(bruneton_ocean freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(bruneton_ocean glut GL)
ENDIF(UNIX AND NOT APPLE)

target_link_libraries(bruneton_ocean glew)
target_link_libraries(bruneton_ocean AntTweakBar)

target_compile_definitions(bruneton_ocean PRIVATE SAMPLE_NAME="Bruneton-Ocean")
