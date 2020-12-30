Multi-Resolution Screen-Space Ambient Occlusion
Thai-Duong Hoang and Kok-Lim Low

You are welcome to use and modify this code for any purpose, but please note that it comes without any guarantees.
Please keep in mind that this project should only be considered a prototype to prove that an idea works. As such, do not expect high code quality.

This project is supposed to be compiled in Release mode using Visual Studio 2010. You can compile and run it directly from the IDE as the paths to shader and mesh files are hardcoded inside main.cpp. If you wish to do otherwise, please modify the paths accordingly.

Once the program is run, use W A S D and the arrow keys to move/rotate the camera. Press F11 to save screenshots (in \img) and Insert to print the camera's orientation and position to a text file (in \txt).

The input mesh is given as a command line parameter (can be configured in Configuration Properties -> Debugging). Three mesh files are given. You can download some very nice scene models at http://www.3drender.com/challenges/index.htm. In main.cpp I've provided some camera and scale parameters for displaying some of those scenes.

You can experiment with the following parameters (found at the beginning of main.cpp, sorry for not providing an easy interface for this):
  - LEVEL_COUNT (default: 5)
  - RESOLUTION (default: 800)
  - scale (default: 1.0) (this should be different for each model)
  - dMax (default: 3.0)
  - rMax (default: 5.0)