#ifndef _RECTANGLE3F_H
#define _RECTANGLE3F_H

#include "vector3f.h"

class Rectangle3f
{
public:
    Rectangle3f();
    Rectangle3f(const Vector3f& v1, const Vector3f& v2);
    Rectangle3f(const Rectangle3f& orig);
    Rectangle3f& operator = (const Rectangle3f& rhs);
    ~Rectangle3f();

    const Vector3f& getV1() const;
    const Vector3f& getV2() const;

    void scale(const Vector3f& v);
    void scale(float s);
    float volume() const;

    void dump() const;

private:
    Vector3f _v1;
    Vector3f _v2;
};

#endif /* end of include guard: rectangle3f.h */