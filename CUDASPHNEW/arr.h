//-----------------------------------------------------------------------------
//  arr.h
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------
//  Utilites for linear arrays
//-----------------------------------------------------------------------------
#include <iostream>


namespace ARR
{

//! Dumps elements of a linear array to stdout.
template <typename T>
inline void Dump (const T* arr, unsigned int numElements, unsigned int pauseAfter = 0)
{
    for (unsigned int i = 0; i < numElements; i++)
    {
        if (pauseAfter != 0)
        {
            if (i % pauseAfter == 0)
            {
                system("pause");
            }
        }
        
        std::cout << arr[i] << std::endl;   
    }
}

//! checks if all elements in a linear array are equal to zero.
template <typename T>
inline bool CheckAllZero (const T* arr, unsigned int numElements)
{
    for (unsigned int i = 0; i < numElements; i++)
    {
        if (arr[i] != static_cast<T>(0))
        {
            return false;
        }
    }

    return true;
}


template <typename T>
inline unsigned int CountElement (const T* arr, unsigned int numElements, 
    const T& element)
{
    unsigned int count = 0;

    for (unsigned int i = 0; i < numElements; i++)
    {
        if (arr[i] == element)
        {
            count++;
        }
    }

    return count;
}

}