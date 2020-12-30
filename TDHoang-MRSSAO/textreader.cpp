#include "textreader.h"

// Prints an error message and terminates the program
void PrintError(const char* message)
{
  printf("%s\n", message);
  //exit(EXIT_FAILURE);
}

// Reads the content of a text file into a null-terminated string
char* ReadFile(const char* file_name)
{
  FILE *fp = fopen(file_name, "r");

  if (fp == NULL)
    PrintError("Cannot open file");

  fseek(fp, 0, SEEK_END);
  int count = ftell(fp);

  char *buffer = (char*) malloc(sizeof(char) * (count + 1));

  rewind(fp);
  count = fread(buffer, sizeof(char), count, fp);
  buffer[count] = '\0';

  fclose(fp);

  return buffer;
}