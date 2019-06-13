SET(SRC_FAST_MASS_SPRING
	./source/timer_wrapper.h
	./source/timer.h
	./source/stb_image_write.h
	./source/stb_image_write.c
	./source/stb_image.h
	./source/stb_image.c
	./source/simulation.h
	./source/simulation.cpp
	./source/scene.h
	./source/scene.cpp
	./source/primitive.h
	./source/primitive.cpp
	./source/opengl_headers.h
	./source/mesh.h
	./source/mesh.cpp
	./source/math_headers.h
	./source/math_headers.cpp
	./source/main.cpp
	./source/io_mesh.h
	./source/io_mesh.cpp
	./source/glsl_wrapper.h
	./source/glsl_wrapper.cpp
	./source/global_headers.h
	./source/fps.h
	./source/fps.cpp
	./source/constraint.h
	./source/constraint.cpp
	./source/clock.h
	./source/camera.h
	./source/camera.cpp
	./source/anttweakbar_wrapper.h
	./source/anttweakbar_wrapper.cpp
)

add_executable(fastmassspring ${SRC_FAST_MASS_SPRING})
target_include_directories(fastmassspring PRIVATE ${PROJECT_SOURCE_DIR}/_libs/eigen)
target_include_directories(fastmassspring PRIVATE ${PROJECT_SOURCE_DIR}/_libs/eigen/Eigen)

IF(WIN32)
target_link_libraries(fastmassspring freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(fastmassspring glut GL)
ENDIF(UNIX AND NOT APPLE)

target_link_libraries(fastmassspring glew)
target_link_libraries(fastmassspring AntTweakBar)
target_link_libraries(fastmassspring tinyxml)

target_compile_definitions(fastmassspring PRIVATE SAMPLE_NAME="FastMassSpring")
