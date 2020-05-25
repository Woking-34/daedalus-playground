set(CUDA_PROJECTIVE_ARMADILLO_SRC
	./CUDA_Projective_Armadillo/ARMADILLO.h
	./CUDA_Projective_Armadillo/OPENGL_DRIVER.h
	./CUDA_Projective_Armadillo/kernel.cu
)

CUDA_ADD_EXECUTABLE(xian_galerkinmg_octopus ${CUDA_PROJECTIVE_ARMADILLO_SRC})
target_link_libraries(xian_galerkinmg_octopus freeglut)
target_link_libraries(xian_galerkinmg_octopus glew)
target_link_libraries(xian_galerkinmg_octopus ${CUDA_LIBRARIES})
target_link_libraries(xian_galerkinmg_octopus ${CUDA_cublas_LIBRARY})
target_link_libraries(xian_galerkinmg_octopus ${CUDA_cusparse_LIBRARY})
target_link_libraries(xian_galerkinmg_octopus ${CUDA_cusolver_LIBRARY})

target_compile_definitions(xian_galerkinmg_octopus PRIVATE SAMPLE_NAME="Xian-GalerkinMG")