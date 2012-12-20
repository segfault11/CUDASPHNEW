#include "vector3f.h"
#include <iostream>
#include <cmath>

Vector3f::Vector3f(): _x(0.0f), _y(0.0f), _z(0.0f) {}
Vector3f::Vector3f(float x, float y, float z): _x(x), _y(y), _z(z) {}
Vector3f::Vector3f(const Vector3f& orig): _x(orig._x), _y(orig._y), 
    _z(orig._z) {}
Vector3f::~Vector3f() {}

Vector3f& Vector3f::operator = (const Vector3f& rhs)
{
    this->_x = rhs._x;
    this->_y = rhs._y;
    this->_z = rhs._z;

    return *this;
}

float Vector3f::getX() const
{
    return _x;
}

float Vector3f::getY() const
{
    return _y;
}

float Vector3f::getZ() const
{
    return _z;
}

void Vector3f::scale(float s)
{
    _x *= s;
    _y *= s;
    _z *= s;
}

void Vector3f::addOffset(float s)
{
    _x += s;
    _y += s;
    _z += s;
}

void Vector3f::add(const Vector3f& v)
{
    this->_x += v._x;
    this->_y += v._y;
    this->_z += v._z;
}

void Vector3f::ceil()
{
    _x = std::ceil(_x);
    _y = std::ceil(_y);
    _z = std::ceil(_z);
}

void Vector3f::dump() const
{
    using std::cout;
    using std::endl;

    cout << _x << endl;
    cout << _y << endl;
    cout << _z << endl;
}

void Vector3f::plus(const Vector3f& v0, const Vector3f& v1, 
        Vector3f& result)
{
    result._x = v0._x + v1._x;
    result._y = v0._y + v1._y;
    result._z = v0._z + v1._z;
}

void Vector3f::minus(const Vector3f& v0, const Vector3f& v1, 
        Vector3f& result)
{
    result._x = v0._x - v1._x;
    result._y = v0._y - v1._y;
    result._z = v0._z - v1._z;
}