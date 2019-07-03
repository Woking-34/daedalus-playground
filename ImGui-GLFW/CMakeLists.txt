add_executable(imgui_glfw
	main.cpp
	
	imgui_impl_glfw.h
	imgui_impl_glfw.cpp
	imgui_impl_opengl2.h
	imgui_impl_opengl2.cpp
)

find_package(OpenGL 4)

target_include_directories(imgui_glfw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries(imgui_glfw ${OPENGL_LIBRARIES})
target_link_libraries(imgui_glfw glfw)
target_link_libraries(imgui_glfw glew)
target_link_libraries(imgui_glfw imgui)