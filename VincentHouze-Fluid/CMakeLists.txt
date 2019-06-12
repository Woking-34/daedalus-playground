file(GLOB_RECURSE SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/Utils/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/CudaCommon/CudaParticlesSystem/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/CudaCommon/CudaParticlesSystem/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/CudaCommon/CudaFluidSolver3D/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/CudaCommon/CudaFluidSolver3D/*.cpp
	
	timer.h
	timer.cpp
	
	nvMath.h
	nvMatrix.h
	nvVector.h
	nvQuaternion.h
	cutil_math.h
)

FILE(GLOB_RECURSE CUFILES "*.cu")

CUDA_ADD_EXECUTABLE(vincenthouze_fluid main.cpp ${SOURCES} ${CUFILES})
target_link_libraries(vincenthouze_fluid freeglut)
target_link_libraries(vincenthouze_fluid glew)
target_link_libraries(vincenthouze_fluid ${CUDA_LIBRARIES})
target_include_directories(vincenthouze_fluid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})