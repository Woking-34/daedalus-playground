CXXFLAGS=-g -Wall -DNDEBUG -fopenmp -ffast-math
CXX=g++

all: fsolver
clean:
	rm -f *.o fsolver

display.o: display.cpp fsolver.h
fsolver.o: fsolver.cpp fsolver.h

fsolver: display.o fsolver.o
	$(CXX) -fopenmp display.o fsolver.o -o fsolver -lglut -lpng -lGLU
