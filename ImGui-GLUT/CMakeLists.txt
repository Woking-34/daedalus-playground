add_executable(imgui_glut
	main.cpp
	
	imgui_impl_glut.h
	imgui_impl_glut.cpp
	imgui_impl_opengl2.h
	imgui_impl_opengl2.cpp
)

target_include_directories(imgui_glut PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries(imgui_glut freeglut)
target_link_libraries(imgui_glut glew)
target_link_libraries(imgui_glut imgui)