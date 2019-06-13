#include <GL/glut.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
//#include <png.h>
#include "fsolver.h"

FluidSolver *fsolver = NULL;
int displayVelocities = 0;
bool running = false;
bool button1Down = false;
bool button2Down = false;
int lastX = 0, lastY = 0;
double angleX = 0.0f, angleY = 0.0f;
double distZ = -140;
int displayMode = 0;
int width = 768;
int height = 768;
int frameCounter = 0;

//void screenshot(int i);
//void exportDensity(int i);

void renderFunc() {
	glDepthMask(GL_TRUE);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	int gridSizeX = fsolver->getGridSizeX(),
	    gridSizeY = fsolver->getGridSizeY(),
		gridSizeZ = fsolver->getGridSizeZ();
	Float dx = fsolver->getVoxelSize();
	const vector<Density> &density = fsolver->getDensity();
	const Vector &curlMagnitude = fsolver->getCurlMagnitude();
	const bool *solid = fsolver->getSolid();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,-distZ,0,0,0,0,-1,0);
	glRotated(angleY, 1, 0, 0);
	glRotated(angleX, 0, 1, 0);
	glTranslatef(-gridSizeX/2.0f, -gridSizeY/2.0f, -gridSizeZ/2.0f);
	glCullFace(GL_NONE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glBegin(GL_QUADS);
		glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		glVertex3i(0,   0, 0);
		glVertex3i(0,   gridSizeY, 0);
		glVertex3i(gridSizeX, gridSizeY, 0);
		glVertex3i(gridSizeX, 0, 0);
		
		glVertex3i(0,   0, gridSizeZ);
		glVertex3i(0,   gridSizeY, gridSizeZ);
		glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
		glVertex3i(gridSizeX, 0, gridSizeZ);
		
		glVertex3i(0, 0,   0);
		glVertex3i(0, 0,   gridSizeZ);
		glVertex3i(0, gridSizeY, gridSizeZ);
		glVertex3i(0, gridSizeY, 0);
		
		glVertex3i(gridSizeX, 0,   0);
		glVertex3i(gridSizeX, 0,   gridSizeZ);
		glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
		glVertex3i(gridSizeX, gridSizeY, 0);

		glVertex3i(0,   0, 0);
		glVertex3i(0,   0, gridSizeZ);
		glVertex3i(gridSizeX, 0, gridSizeZ);
		glVertex3i(gridSizeX, 0, 0);
		
		glVertex3i(0,   gridSizeY, 0);
		glVertex3i(0,   gridSizeY, gridSizeZ);
		glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
		glVertex3i(gridSizeX, gridSizeY, 0);
	glEnd();			
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);
	glColor4f(.5f, .5f, 1.0f, 1.0f);
	for (int z=0, pos = 0; z<gridSizeZ; ++z) {
		for (int y=0; y<gridSizeY; ++y) {
			for (int x=0; x<gridSizeX; ++x, ++pos) {
				if (!solid[pos])
					continue;
				glVertex3i(x,   y, z);
				glVertex3i(x,   y+1, z);
				glVertex3i(x+1, y+1, z);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y, z+1);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y, z+1);
				
				glVertex3i(x, y,   z);
				glVertex3i(x, y,   z+1);
				glVertex3i(x, y+1, z+1);
				glVertex3i(x, y+1, z);
				
				glVertex3i(x+1, y,   z);
				glVertex3i(x+1, y,   z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);

				glVertex3i(x,   y, z);
				glVertex3i(x,   y, z+1);
				glVertex3i(x+1, y, z+1);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y+1, z);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
			}
		}
	}
	glEnd();
	glDepthMask(GL_FALSE);
	glBegin(GL_QUADS);
	for (int z=0, pos = 0; z<gridSizeZ; ++z) {
		for (int y=0; y<gridSizeY; ++y) {
			for (int x=0; x<gridSizeX; ++x, ++pos) {
				if (solid[pos])
					continue;
				if (displayMode == 0) {
					Float f = .5f * density[pos].density;
					glColor4f(1.0f, 1.0f, 1.0f, f);
				} else {
					Float f = .01f * curlMagnitude[pos];
					glColor4f(0.0f, 0.0f, 1.0f, f);
				}
				glVertex3i(x,   y, z);
				glVertex3i(x,   y+1, z);
				glVertex3i(x+1, y+1, z);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y, z+1);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y, z+1);
				
				glVertex3i(x, y,   z);
				glVertex3i(x, y,   z+1);
				glVertex3i(x, y+1, z+1);
				glVertex3i(x, y+1, z);
				
				glVertex3i(x+1, y,   z);
				glVertex3i(x+1, y,   z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);

				glVertex3i(x,   y, z);
				glVertex3i(x,   y, z+1);
				glVertex3i(x+1, y, z+1);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y+1, z);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
			}
		}
	}

	glEnd();
	if (displayVelocities > 0) {
		glDisable(GL_BLEND);
		glBegin(GL_LINES);
		glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
		for (int z=0, pos=0; z<gridSizeZ; ++z) {
			for (int y=0; y<gridSizeY; ++y) {
				for (int x=0; x<gridSizeX; ++x, ++pos) {
					if ((x%2) != 0 || (y%2) != 0 || (z%2) != 0)
						continue;
					Vector3 velocity;
					if (displayVelocities == 1)
						velocity = fsolver->getVelocity(Point3((x+.5f)*dx, (y+.5f)*dx, (z+.5f)*dx));
					else if (displayVelocities == 2)
						velocity = fsolver->getForce(Point3((x+.5f)*dx, (y+.5f)*dx, (z+.5f)*dx));
					velocity /= dx;
					glVertex3f(x+.5f, y+.5f, z+.5f);
					glVertex3f(x+.5f+velocity.x, y+.5f+velocity.y, z+.5f+velocity.z);
				}
			}
		}
		glEnd();
	}

	glutSwapBuffers();
	if (running) {
		fsolver->step(0.1f);
		//exportDensity(frameCounter);
		//screenshot(frameCounter++);
		glutPostRedisplay();
	}
}

//void exportDensity(int i) {
//	char fname[128], cmd[128];
//	snprintf(fname, sizeof(fname), "output/density-%04i.vol", i);
//	snprintf(cmd, sizeof(fname), "gzip -f output/density-%04i.vol&", i);
//	cout << "    + Saving \"" << fname << "\"" << endl;
//	std::ofstream os(fname);
//
//	int xres = fsolver->getGridSizeX(),
//	    yres = fsolver->getGridSizeY(),
//		zres = fsolver->getGridSizeZ();
//	const vector<Density> &density = fsolver->getDensity();
//
//	Float scale = 1.0f / std::max(std::max(xres, yres), zres);
//
//	os.write("VOL", 3);
//	char version = 3;
//	os.write((char *) &version, sizeof(char));
//	int value = 1;
//	os.write((char *) &value, sizeof(int));
//	os.write((char *) &xres, sizeof(int));
//	os.write((char *) &yres, sizeof(int));
//	os.write((char *) &zres, sizeof(int));
//	value = 1;
//	os.write((char *) &value, sizeof(int));
//	
//	float minX = -xres/2.0f*scale;
//	float minY = -yres/2.0f*scale;
//	float minZ = -zres/2.0f*scale;
//	float maxX = xres/2.0f*scale;
//	float maxY = yres/2.0f*scale;
//	float maxZ = zres/2.0f*scale;
//
//	os.write((char *) &minX, sizeof(float));
//	os.write((char *) &minY, sizeof(float));
//	os.write((char *) &minZ, sizeof(float));
//	os.write((char *) &maxX, sizeof(float));
//	os.write((char *) &maxY, sizeof(float));
//	os.write((char *) &maxZ, sizeof(float));
//	for (size_t i=0; i<density.size(); ++i) {
//		float value = density[i].density;
//		os.write((char *) &value, sizeof(float));
//	}
//	os.close();
//	
//	system(cmd);
//}

//void screenshot(int i) {
//	char fname[128];
//	snprintf(fname, sizeof(fname), "output/screenshot-%04i.png", i);
//	cout << "    + Saving \"" << fname << "\"" << endl;
//
//	FILE *file = fopen(fname, "wb");
//	if (file == NULL) {
//		cerr << "Error: Could not open \"" << fname << "\"!" << endl;
//		exit(-1);
//	}
//
//	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//	if (png_ptr == NULL) {
//		cout << "Error: could not initialize libpng!" << endl;
//		exit(-1);
//	}
//
//	png_infop info = png_create_info_struct(png_ptr);
//	if (info == NULL) {
//		cout << "Error: could not initialize libpng!" << endl;
//		exit(-1);
//	}
//
//	if (setjmp(png_ptr->jmpbuf)) {
//		cerr << "Error: \"" << fname << "\" could not be written!" 
//			 << endl;
//		exit(-1);
//		
//	}
//
//	png_init_io (png_ptr, file);
//	png_set_compression_level(png_ptr, 3);
//	png_set_IHDR(png_ptr, info, width, height, 8, PNG_COLOR_TYPE_RGB,
//		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
//	png_write_info(png_ptr, info);
//
//	png_bytepp rows = new png_bytep[height];
//
//	png_bytep data = new png_byte[3*width*height];
//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
//	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);
//	for (int i=0; i<height; i++)
//		rows[height-1-i] = data + width*3*i;
//
//	png_write_image(png_ptr, (png_bytepp) rows);
//	png_write_end(png_ptr, info);
//	png_destroy_write_struct(&png_ptr, &info);
//	
//	delete[] data;
//	delete[] rows;
//	fclose(file);
//}

void mouseFunc(int button, int state, int x, int y) {
	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON) {
		button1Down = true;
		lastX = x;
		lastY = y;
	} else if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) {
		button1Down = false;
	} else if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON) {
		button2Down = true;
		lastX = x;
		lastY = y;
	} else if (state == GLUT_UP && button == GLUT_RIGHT_BUTTON) {
		button2Down = false;
	}
}

void motionFunc(int x, int y) {
	if (button1Down) {
		double dx = (x - lastX)/2.0f;
		double dy = (y - lastY)/2.0f;
		angleX -= dx;
		angleY += dy;
		lastX = x;
		lastY = y;
		glutPostRedisplay();
	} else if (button2Down) {
		double dy = (y - lastY)/10.0f;
		distZ -= dy;
		lastY = y;
		glutPostRedisplay();
	}
}

void reshapeFunc(int w, int h) {
	if (h == 0 || w == 0)
		return;
	glMatrixMode(GL_PROJECTION);
	glViewport(0, 0, w, h);
	glLoadIdentity();
	gluPerspective(45, (float)w/(float) h, 0.1, 1e5);
	glMatrixMode(GL_MODELVIEW);
}

void keyboardFunc(unsigned char key, int x, int y) {
	if (key == 'v')
		displayVelocities = (displayVelocities + 1)%3;
	if (key == ' ')
		displayMode = (displayMode + 1)%2;
	if (key == 's')
		fsolver->step(1.0f);
	if (key == 'r')
		running = !running;
	if (key == 27)
		exit(0);

	glutPostRedisplay();
}

int main(int argc, char **argv) {
	int resX = 80, resY = 80, resZ = 40;
	Float size = 0.5f/60 * 5;
	fsolver = new FluidSolver(resX, resY, resZ, size);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(width,height);
	glutCreateWindow("Fluid Solver [density]");
	glutDisplayFunc(renderFunc);
	glutMouseFunc(mouseFunc);
	glutReshapeFunc(reshapeFunc);
	glutKeyboardFunc(keyboardFunc);
	glutMotionFunc(motionFunc);
	glutMainLoop();
	delete fsolver;
	return 0;
}
