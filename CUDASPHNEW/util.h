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
#include <iostream>

#define M_PI       3.14159265358979323846

#define CUDA_SAFE_CALL(x) cudaSafeCall(x, __FILE__, __LINE__)
#define THROW_EXCEPTION(x) throwException(x, __FILE__, __LINE__)

inline void cudaSafeCall (cudaError_t err, const char* filename, unsigned int line)
{
    using namespace std;

    if (err == cudaSuccess)
    {
        return;
    }    
    
    stringstream str(stringstream::in | stringstream::out);

    str << "Exception thrown in FILE: " << filename << " LINE: " << line << endl;
    str << "Error Message: " << cudaGetErrorString(cudaGetLastError()) << endl;
    
    throw runtime_error(str.str());
}

inline void throwException (const char* message, const char* filename, 
    unsigned int line)
{
    using namespace std;

    stringstream str(stringstream::in | stringstream::out);

    str << "Exception thrown in FILE: " << filename << " LINE: " << line << endl;
    str << "Error Message: " << message << endl;

    throw runtime_error(str.str());
}

template<typename T>
inline void saveDelete (T** ptr)
{
    if (*ptr != NULL) {
        delete *ptr;
        *ptr = NULL;
    }
}

template<typename T>
inline void saveDeleteArray (T** ptr)
{
    if (*ptr != NULL) {
        delete *ptr;
        *ptr = NULL;
    }
}

template<typename T>
void cudaSafeFree (T** ptr ) 
{
    if (*ptr != NULL) {
        cudaFree(*ptr);
        *ptr = NULL;
    }
}


template<typename T> 
inline void CUDADumpArray (T* arr, unsigned int numElements, unsigned int offset = 0, 
    unsigned int stride = 1, unsigned int pauseAfter = 0)
{
    T* hostData = new T[numElements];

    CUDA_SAFE_CALL( cudaMemcpy(hostData, arr, sizeof(T)*numElements,
        cudaMemcpyDeviceToHost) );
    
    for (unsigned int i = 0; i < numElements; i++)
    {
        std::cout << hostData[i*stride + offset] << std::endl;

        if (pauseAfter != 0)
        {
            if ((i % pauseAfter) == 0)
            {
                system("pause");
            }
        }
    }
    
    delete[] hostData;
}


template<typename T> 
inline void CUDADumpArrayRanged (T* arr, unsigned int numElements, 
    unsigned int s, unsigned int e, unsigned int offset = 0, 
    unsigned int stride = 1, unsigned int pauseAfter = 0)
{
    T* hostData = new T[numElements];

    CUDA_SAFE_CALL( cudaMemcpy(hostData, arr, sizeof(T)*numElements,
        cudaMemcpyDeviceToHost) );
    
    for (unsigned int i = s; i <= e; i++)
    {
        std::cout << hostData[i*stride + offset] << std::endl;

        if (pauseAfter != 0)
        {
            if ((i % pauseAfter) == 0)
            {
                system("pause");
            }
        }
    }
    
    delete[] hostData;
}

template<typename T>
void CUDADumpArrayElements (T* arr, unsigned int numElemsArr, 
    int* elemArr, unsigned int numElements, unsigned int offset = 0, 
    unsigned int stride = 1, unsigned int pauseAfter = 0)
{
    T* hostData = new T[numElemsArr];
    int* elemData = new int[numElements];

    CUDA_SAFE_CALL( cudaMemcpy(hostData, arr, sizeof(T)*numElemsArr,
        cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(elemData, elemArr, sizeof(int)*numElements,
        cudaMemcpyDeviceToHost) );
    
    for (unsigned int i = 0; i < numElements; i++)
    {
        std::cout << hostData[elemData[i]*stride + offset] << std::endl;
        
        if (pauseAfter != 0)
        {
            if ((i % pauseAfter) == 0)
            {
                system("pause");
            }
        } 
    }
    
    delete[] hostData;
    delete[] elemData;
}

inline void cudaCheckUnique (unsigned int* arr, unsigned int numElements, unsigned int maxRange)
{
    bool* test = new bool[maxRange + 1];
    int* hostArr = new int[numElements];
    
    for (unsigned int i = 0; i < maxRange + 1; i++)
    {
        test[i] = false;
    }

    CUDA_SAFE_CALL( cudaMemcpy(hostArr, arr, numElements*sizeof(int), cudaMemcpyDeviceToHost) ); 
    
    for (unsigned int i = 0; i < numElements; i++)
    {
        if (test[hostArr[i]] == true)
        {
            std::cout << "at least one element exists at least twice in the array." << std::endl;
            system("pause");
        }

        test[hostArr[i]] = true;
    }

    std::cout << "all elements are unique" << std::endl;

    delete[] hostArr;
    delete[] test;
}

class CudaTimer
{
    enum 
    {
        CT_TAKING_TIME,
        CT_STOPPED
    };

public:
    CudaTimer();
    ~CudaTimer();
    void Start();
    void Stop();
    void DumpElapsed() const;

private:
    cudaEvent_t mStart;
    cudaEvent_t mStop;
    float mTime;
    int mState;
};

#endif /* end of include guard: util.h */