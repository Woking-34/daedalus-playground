#ifndef __textreader_h__
#define __textreader_h__

#include <stdio.h>
#include <stdlib.h>

// Prints an error message and terminates the program
void PrintError(const char* message);

// Reads the content of a text file into a null-terminated string
char* ReadFile(const char* file_name);

#endif __textreader_h__