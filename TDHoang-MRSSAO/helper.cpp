#include "helper.h"

void PrintGLError()
{
  GLenum errorCode = glGetError();
  while (errorCode != GL_NO_ERROR)
  {
    printf("OpenGL error: %s\n", gluErrorString(errorCode));
    errorCode = glGetError();
  }
}