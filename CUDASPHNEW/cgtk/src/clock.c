#include <Windows.h>
#include <stdio.h>
#include "../include/clock.h"

static LARGE_INTEGER gStart;


void cgtkClockStart() {
	QueryPerformanceCounter(&gStart);
}

double cgtkClockElapsed() {
	LARGE_INTEGER frequency;     
	LARGE_INTEGER end;           
   
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&end);

	return (end.QuadPart - gStart.QuadPart)*1000.0f/frequency.QuadPart;
}

void cgtkClockDumpElapsed() {
	printf("elapsed: %lf ms\n", cgtkClockElapsed());
}