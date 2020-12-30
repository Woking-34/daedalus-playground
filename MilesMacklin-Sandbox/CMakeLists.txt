set(MMACKLIN_CORE_SRC
	core/aabbtree.cpp
	core/aabbtree.h
	core/core.cpp
	core/core.h
	core/hashgrid.cpp
	core/hashgrid.h
	core/log.cpp
	core/log.h
	core/mat22.h
	core/mat33.h
	core/mat44.h
	core/maths.cpp
	core/maths.h
	core/memory.cpp
	core/memory.h
	core/mesh.cpp
	core/mesh.h
	core/meshuv.cpp
	core/meshuv.h
	core/perlin.cpp
	core/perlin.h
	core/pfm.cpp
	core/pfm.h
	core/platform.cpp
	core/platform.h
	core/png.cpp
	core/png.h
	core/point3.h
#	core/radix.cpp
#	core/radix.h
	core/sdf.cpp
	core/sdf.h
	core/shader.cpp
	core/shader.h
	core/skylight.cpp
	core/skylight.h
	core/tga.cpp
	core/tga.h
	core/threadgroup.cpp
	core/threadgroup.h
	core/threading.cpp
	core/threading.h
#	core/threadinglinux.h
#	core/threadingwin32.h
	core/types.h
	core/vec2.h
	core/vec3.h
	core/vec4.h
	core/voxelize.cpp
	core/voxelize.h
)

add_library(mmacklin_core STATIC ${MMACKLIN_CORE_SRC})
target_include_directories(mmacklin_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mmacklin_sh projects/sh/main.cpp projects/sh/sh.cpp projects/sh/sh.h)
target_link_libraries(mmacklin_sh mmacklin_core)
target_link_libraries(mmacklin_sh glew)
target_compile_definitions(mmacklin_sh PRIVATE SAMPLE_NAME="MilesMacklin-Sandbox")

IF(WIN32)
target_link_libraries(mmacklin_sh freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(mmacklin_sh glut GLU GL)
ENDIF(UNIX AND NOT APPLE)

add_executable(mmacklin_smoke projects/smoke/fluid3d.cpp projects/smoke/blackbody.cpp projects/smoke/blackbody.h projects/smoke/grid3d.h)
target_link_libraries(mmacklin_smoke mmacklin_core)
target_link_libraries(mmacklin_smoke glew)
target_compile_definitions(mmacklin_smoke PRIVATE SAMPLE_NAME="MilesMacklin-Sandbox")

IF(WIN32)
target_link_libraries(mmacklin_smoke freeglut)
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
target_link_libraries(mmacklin_smoke glut GLU GL)
ENDIF(UNIX AND NOT APPLE)
