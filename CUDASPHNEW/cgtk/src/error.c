#include "../include/error.h"
#include <stdio.h>

void cgtkDumpErrorMessage(const char* filename, unsigned int linenumber,
    const char* message)
{
	printf("Error in FILE: %s, LINE: %d\nMessage: %s\n", filename, linenumber,
		message);
}

