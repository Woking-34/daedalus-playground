SET(SRC_BULLET_PRIMITIVESHAPES
	bullet_PrimitiveShapes.cpp
	
	geometry.h
)

add_executable( bullet_primitiveshapes ${SRC_BULLET_PRIMITIVESHAPES} )

IF(WIN32)
	target_link_libraries( bullet_primitiveshapes glew freeglut ${BULLET_LIBRARIES} )
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
	target_link_libraries( bullet_primitiveshapes GLEW glut GL ${BULLET_LIBRARIES} )
ENDIF(UNIX AND NOT APPLE)

target_include_directories(bullet_primitiveshapes PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${BULLET_INCLUDE_DIRS})