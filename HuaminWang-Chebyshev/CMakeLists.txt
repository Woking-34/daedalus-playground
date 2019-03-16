#set(HUAMINWANG_CHEBYSHEV_SRC
#	TIMER.h
#	TET_MESH.h
#	PROJECTIVE_TET_MESH.h
#	PROJECTIVE_MESH.h
#	MY_MATH.h
#	MY_GLSL.h
#	MESH.h
#	IO_FUNC.h
#	INTERSECTION.h
#	DYNAMIC_MESH.h
#	DISTANCE.h
#	CUDA_PROJECTIVE_TET_MESH.h
#	CUDA_PROJECTIVE_MESH.h
#	BMP_IO.h
#)

set(CUDA_PROJECTIVE_ARMADILLO_SRC
	./CUDA_Projective_Armadillo/ARMADILLO.h
	./CUDA_Projective_Armadillo/OPENGL_DRIVER.h
	./CUDA_Projective_Armadillo/kernel.cu
)

set(CUDA_PROJECTIVE_HAMMOCK_SRC
	./CUDA_Projective_Hammock/HAMMOCK.h
	./CUDA_Projective_Hammock/OPENGL_DRIVER.h
	./CUDA_Projective_Hammock/kernel.cu
)

set(CUDA_PROJECTIVE_SQUARE_SRC
	./CUDA_Projective_Square/CLOTHING.h
	./CUDA_Projective_Square/OPENGL_DRIVER.h
	./CUDA_Projective_Square/kernel.cu
)

CUDA_ADD_EXECUTABLE(wang_chebyshev_armadillo ${CUDA_PROJECTIVE_ARMADILLO_SRC})
target_link_libraries(wang_chebyshev_armadillo freeglut)
target_link_libraries(wang_chebyshev_armadillo glew)
target_link_libraries(wang_chebyshev_armadillo ${CUDA_LIBRARIES})

CUDA_ADD_EXECUTABLE(wang_chebyshev_hammock ${CUDA_PROJECTIVE_HAMMOCK_SRC})
target_link_libraries(wang_chebyshev_hammock freeglut)
target_link_libraries(wang_chebyshev_hammock glew)
target_link_libraries(wang_chebyshev_hammock ${CUDA_LIBRARIES})

CUDA_ADD_EXECUTABLE(wang_chebyshev_square ${CUDA_PROJECTIVE_SQUARE_SRC})
target_link_libraries(wang_chebyshev_square freeglut)
target_link_libraries(wang_chebyshev_square glew)
target_link_libraries(wang_chebyshev_square ${CUDA_LIBRARIES})

target_compile_definitions(wang_chebyshev_armadillo PRIVATE SAMPLE_NAME="HuaminWang-Chebyshev")
target_compile_definitions(wang_chebyshev_hammock PRIVATE SAMPLE_NAME="HuaminWang-Chebyshev")
target_compile_definitions(wang_chebyshev_square PRIVATE SAMPLE_NAME="HuaminWang-Chebyshev")