#include "util.h"
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------------
CudaTimer::CudaTimer ()
    : mState(CT_STOPPED)
{
    cudaEventCreate(&mStart);
    cudaEventCreate(&mStop);
}
//-----------------------------------------------------------------------------
CudaTimer::~CudaTimer ()
{
    cudaEventDestroy(mStart);
    cudaEventDestroy(mStop);
}
//-----------------------------------------------------------------------------
void CudaTimer::Start () 
{
    cudaEventRecord(mStart, 0);
    mState = CT_TAKING_TIME;
}
//-----------------------------------------------------------------------------
void CudaTimer::Stop ()
{
    cudaEventRecord(mStop, 0);
    cudaEventSynchronize(mStop);
    cudaEventElapsedTime(&mTime, mStart, mStop);
    mState = CT_STOPPED;
}
//-----------------------------------------------------------------------------
void CudaTimer::DumpElapsed () const
{
    if (mState == CT_TAKING_TIME)
    {
        return;
    }

    cout << mTime << " [ms] elapsed" << endl;
}
