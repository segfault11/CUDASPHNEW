#include "../include/math.h"
#include <math.h>

void cgtkVec3fDotProduct(CGTKVec3f u, CGTKVec3f v, float* result) {
	*result = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void cgtkVec3fCrossProduct(CGTKVec3f u, CGTKVec3f v, CGTKVec3f w) {
	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - u[0]*v[2];
	w[2] = u[0]*v[1] - u[1]*v[0];
}

void cgtkVec3fNormalize(CGTKVec3f v) {
	float norm = sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

