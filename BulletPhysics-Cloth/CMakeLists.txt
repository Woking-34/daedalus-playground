SET(SRC_BULLET_CLOTH
	bullet_Cloth.cpp
	
	geometry.h
)

add_executable( bullet_cloth ${SRC_BULLET_CLOTH} )

IF(WIN32)
	target_link_libraries( bullet_cloth glew freeglut ${BULLET_LIBRARIES} )
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
	target_link_libraries( bullet_cloth glew glut GLU GL ${BULLET_LIBRARIES} )
ENDIF(UNIX AND NOT APPLE)

target_include_directories(bullet_cloth PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${BULLET_INCLUDE_DIRS})
