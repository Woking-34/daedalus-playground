set(CUDA_RAYTRACINGGEMS_VOLUMEPT_SRC
	hdr_loader.h
	
	main.cpp
	
	volume_kernel.h
	volume_kernel.cu
)

find_package(OpenGL 4)

CUDA_ADD_EXECUTABLE(raytracinggems_volumept ${CUDA_RAYTRACINGGEMS_VOLUMEPT_SRC})
target_link_libraries(raytracinggems_volumept ${OPENGL_LIBRARIES})
target_link_libraries(raytracinggems_volumept glew)
target_link_libraries(raytracinggems_volumept glfw)
target_link_libraries(raytracinggems_volumept ${CUDA_LIBRARIES})
target_include_directories(raytracinggems_volumept PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})