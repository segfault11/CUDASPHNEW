#ifndef _VECTOR3F_H
#define _VECTOR3F_H

class Rectangle3f;

class Vector3f 
{
    friend Rectangle3f;
public:
    Vector3f();
    Vector3f(float x, float y, float z);
    Vector3f(const Vector3f& orig);
    Vector3f& operator = (const Vector3f& rhs);
    ~Vector3f();

    float getX() const;
    float getY() const;
    float getZ() const;

    void scale(float s);
    void addOffset(float s);

    void add(const Vector3f& v);

    void ceil();

    void dump() const;

    static void minus(const Vector3f& v0, const Vector3f& v1, 
        Vector3f& result);

    static void plus(const Vector3f& v0, const Vector3f& v1, 
        Vector3f& result);
private:
    float _x;
    float _y;
    float _z;
};

#endif /* end of include guard: vector3f.h */ 