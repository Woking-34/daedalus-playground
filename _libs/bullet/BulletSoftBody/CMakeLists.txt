SET(BulletSoftBody_SRCS
	btSoftBody.cpp
	btSoftBodyConcaveCollisionAlgorithm.cpp
	btSoftBodyHelpers.cpp
	btSoftBodyRigidBodyCollisionConfiguration.cpp
	btSoftRigidCollisionAlgorithm.cpp
	btSoftRigidDynamicsWorld.cpp
	btSoftSoftCollisionAlgorithm.cpp
	btDefaultSoftBodySolver.cpp
)

SET(BulletSoftBody_HDRS
	btSoftBody.h
	btSoftBodyData.h
	btSoftBodyConcaveCollisionAlgorithm.h
	btSoftBodyHelpers.h
	btSoftBodyRigidBodyCollisionConfiguration.h
	btSoftRigidCollisionAlgorithm.h
	btSoftRigidDynamicsWorld.h
	btSoftSoftCollisionAlgorithm.h
	btSparseSDF.h
	
	btSoftBodySolvers.h
	btDefaultSoftBodySolver.h
	
	btSoftBodySolverVertexBuffer.h
)

ADD_LIBRARY(BulletSoftBody STATIC ${BulletSoftBody_SRCS} ${BulletSoftBody_HDRS})
TARGET_LINK_LIBRARIES(BulletSoftBody BulletDynamics)

SET_PROPERTY(TARGET BulletSoftBody PROPERTY FOLDER "_library/bullet")
