#ifndef __shader_h__
#define __shader_h__

#include "GL/glew.h"
#include "GL/freeglut.h"
#include "textreader.h"
#include "helper.h"

// Shader types
typedef enum
{
  EVertexShader,
  EFragmentShader,
} EShaderType;

void AttachShader(GLuint program, GLuint shader);
void LinkProgram(GLuint program);
GLuint CreateProgram();

// Loads a shader from a file
GLuint LoadShader(EShaderType type, char* file_name, std::string rootData);

// Prints a shader object's info log
void PrintShaderLog(GLuint shader);

// Prints a program object's info log
void PrintProgramLog(GLuint program);

#endif __shader_h__