cmake_minimum_required(VERSION 2.8.11)
cmake_policy(VERSION 2.8.11)

project(bullet)

SET(BULLET_LIBRARIES
	LinearMath
	BulletCollision
	BulletDynamics
	BulletSoftBody
	
	PARENT_SCOPE
)

SET(BULLET_INCLUDE_DIRS
	${CMAKE_CURRENT_LIST_DIR}
	
	PARENT_SCOPE
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})

ADD_SUBDIRECTORY( LinearMath )
ADD_SUBDIRECTORY( BulletCollision )
ADD_SUBDIRECTORY( BulletDynamics )
ADD_SUBDIRECTORY( BulletSoftBody )