//-----------------------------------------------------------------------------
//  Point3f.h
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------
#ifndef _POINT3F_H
#define _POINT3F_H

class Point3f
{
public:
    Point3f ();
    Point3f (float x, float y, float y);
    Point3f (const Point3f& orig);
    ~Point3f ();

    Point3f& operator= (const Point3f& orig);

    // class access
    inline float GetX () const;
    inline float GetY () const;
    inline float GetZ () const;

    inline void SetX (float x);
    inline void SetY (float y);
    inline void SetZ (float z);


private:
    float mX, mY, mZ;
};

#include "Point3f.inl

#endif // Point3f.h