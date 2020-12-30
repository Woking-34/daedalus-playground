SET(SRC_MRSSAO_WHITECAPS
	frame3d.cpp
	frame3d.h
	frameGrab.cpp
	frameGrab.h
	frameRate.cpp
	frameRate.h
	glm.cpp
	glm.h
	helper.cpp
	helper.h
	HighResolutionTimer.h

	math3d.cpp
	math3d.h
	shader.cpp
	shader.h
	textreader.cpp
	textreader.h
	
	main.cpp
)

add_executable(tdhoang_mrssao ${SRC_MRSSAO_WHITECAPS})

IF(WIN32)
target_link_libraries(tdhoang_mrssao freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(tdhoang_mrssao glut GLU GL)
ENDIF(UNIX AND NOT APPLE)

target_link_libraries(tdhoang_mrssao glew)

target_compile_definitions(tdhoang_mrssao PRIVATE SAMPLE_NAME="TDHoang-MRSSAO")
