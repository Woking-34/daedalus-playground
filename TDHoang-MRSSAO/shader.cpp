#include "shader.h"
#include "helper.h"

GLuint CreateProgram()
{
  return glCreateProgram();
}

void AttachShader(GLuint program, GLuint shader)
{
  glAttachShader(program, shader);
}

void LinkProgram(GLuint program)
{
  glLinkProgram(program);
  int program_linked;
  glGetProgramiv(program, GL_LINK_STATUS, &program_linked);
  if (!program_linked)
    PrintError("Program cannot be linked");

  PrintProgramLog(program);

  //glUseProgram(program);
}

// Loads a shader from a file
GLuint LoadShader(EShaderType type, char *file_name, std::string rootData)
{
  char *temp = ReadFile((rootData + file_name).c_str());
  const char *shader_text = temp;

  GLuint shader;
  if (type == EVertexShader)
    shader = glCreateShader(GL_VERTEX_SHADER);
  else if (type == EFragmentShader)
    shader = glCreateShader(GL_FRAGMENT_SHADER);
  else
    PrintError("Wrong shader type");

  glShaderSource(shader, 1, &shader_text, NULL);

  glCompileShader(shader);
  int shader_compiled;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_compiled);
  if (!shader_compiled)
    PrintError("Shader cannot be compiled");

  PrintShaderLog(shader);

  free(temp);

  return shader;
}

// Prints a shader object's info log
void PrintShaderLog(GLuint shader)
{
  int log_length = 0;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_length);

  if (log_length > 0)
  {
    char *log = (char*) malloc(sizeof(char) * log_length);
    glGetShaderInfoLog(shader, log_length, NULL, log);
    printf("Shader log:\n");
    printf("%s\n", log);
    free(log);
  }
}

// Prints a program object's info log
void PrintProgramLog(GLuint program)
{
  int log_length = 0;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &log_length);

  if (log_length > 0)
  {
    char *log = (char*) malloc(sizeof(char) * log_length);
    glGetProgramInfoLog(program, log_length, NULL, log);
    printf("Program log:\n");
    printf("%s\n", log);
    free(log);
  }
}