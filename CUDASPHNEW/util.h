#ifndef UTIL_H
#define UTIL_H

#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include <gl\freeglut_ext.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <exception>
#include <stdio.h>      
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>

#define M_PI       3.14159265358979323846

#define CUDA_SAFE_CALL(x) cudaSafeCall(x, __FILE__, __LINE__)
#define THROW_EXCEPTION(x) throwException(x, __FILE__, __LINE__)

inline void cudaSafeCall(cudaError_t err, const char* filename, unsigned int line)
{
    using namespace std;

    if (err == cudaSuccess) {
        return;
    }    
    
    stringstream str(stringstream::in | stringstream::out);

    str << "Exception thrown in FILE: " << filename << " LINE: " << line << endl;
    str << "Error Message: " << cudaGetErrorString(cudaGetLastError()) << endl;
    
    throw runtime_error(str.str());
}

inline void throwException(const char* message, const char* filename, 
    unsigned int line)
{
    using namespace std;

    stringstream str(stringstream::in | stringstream::out);

    str << "Exception thrown in FILE: " << filename << "LINE: " << line << endl;
    str << "Error Message: " << message << endl;

    throw runtime_error(str.str());
}

template<typename T>
inline void saveDelete(T** ptr)
{
    if (*ptr != NULL) {
        delete *ptr;
        *ptr = NULL;
    }
}

template<typename T>
inline void saveDeleteArray(T** ptr)
{
    if (*ptr != NULL) {
        delete *ptr;
        *ptr = NULL;
    }
}

template<typename T>
void cudaSafeFree(T** ptr ) 
{
    if (*ptr != NULL) {
        cudaFree(*ptr);
        *ptr = NULL;
    }
}


#endif /* end of include guard: util.h */