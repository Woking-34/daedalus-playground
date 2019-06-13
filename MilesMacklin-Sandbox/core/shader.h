#pragma once

#include <GL/glew.h>
#include <GL/freeglut.h>

#include "core/maths.h"

#define glVerify(x) {x; glAssert(#x, __LINE__, __FILE__);}
void glAssert(const char* msg, long line, const char* file);

GLuint CompileProgramFromFile(const char *vertexPath, const char *fragmentPath);
GLuint CompileProgram(const char *vsource=NULL, const char *fsource=NULL, const char* gsource=NULL);
GLuint CompileProgram(const char *vsource, const char* csource, const char* esource, const char* fsource);

void DrawPlane(const Vec4& p, bool color=true);
void DrawString(int x, int y, const char* s, ...);
void DrawFrustum(const Matrix44& projToWorld);
