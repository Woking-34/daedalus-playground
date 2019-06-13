SET(SRC_WENZELJAKOB_FLUID
	aabb.h
	common.h
	sphere.h
	vector.h
	
	fsolver.h
	fsolver.cpp
	display.cpp
)

find_package(Boost)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
	
	add_executable(wenzeljakob_fluid ${SRC_WENZELJAKOB_FLUID})
	
	IF(WIN32)
	target_link_libraries(wenzeljakob_fluid freeglut)
	ENDIF(WIN32)
	IF(UNIX AND NOT APPLE)
	target_link_libraries(wenzeljakob_fluid glut GL)
	ENDIF(UNIX AND NOT APPLE)
	
	IF(WIN32)
		set_target_properties(wenzeljakob_fluid PROPERTIES COMPILE_FLAGS "/openmp")
	ENDIF(WIN32)
	IF(UNIX AND NOT APPLE)
		set_target_properties(wenzeljakob_fluid PROPERTIES COMPILE_FLAGS "-fopenmp")
		set_target_properties(wenzeljakob_fluid PROPERTIES LINK_FLAGS "-fopenmp") 
	ENDIF(UNIX AND NOT APPLE)
	
	target_compile_definitions(wenzeljakob_fluid PRIVATE SAMPLE_NAME="WenzelJakob-Fluid")
endif()