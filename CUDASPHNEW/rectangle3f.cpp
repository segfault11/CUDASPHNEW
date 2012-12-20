#include <iostream>
#include "rectangle3f.h"
#include <cmath>

Rectangle3f::Rectangle3f() {}
Rectangle3f::Rectangle3f(const Vector3f& v1, const Vector3f& v2):
    _v1(v1), _v2(v2) {}
Rectangle3f::Rectangle3f(const Rectangle3f& orig): _v1(orig._v1), 
    _v2(orig._v2) {};
Rectangle3f::~Rectangle3f() {}

Rectangle3f& Rectangle3f::operator = (const Rectangle3f& rhs) 
{
    this->_v1 = rhs._v1;
    this->_v2 = rhs._v2;

    return *this;
}

const Vector3f& Rectangle3f::getV1() const
{
    return _v1;
}

const Vector3f& Rectangle3f::getV2() const
{
    return _v2;
}

void Rectangle3f::scale(float s)
{
    _v1.scale(s);
    _v2.scale(s);
}

float Rectangle3f::volume() const
{
    return std::abs((_v1.getX() - _v2.getX())*(_v1.getY() - _v2.getY())
        *(_v1.getZ() - _v2.getZ()));
}

void Rectangle3f::dump() const
{
    using std::cout;
    using std::endl;

    cout << "v1" << endl;
    _v1.dump();
    cout << "v2" << endl;
    _v2.dump();
}