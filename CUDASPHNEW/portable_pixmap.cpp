#include "portable_pixmap.h"
#include <fstream>
#include <cmath>
#include <fstream>

inline float clamp(float x, float min, float max);

PortablePixmap::PortablePixmap(unsigned int width, unsigned int height, 
    unsigned int maxColor) : _width(width), _height(height), 
    _maxColor(maxColor)
{
    _data = new unsigned int[3*_width*_height];
    memset(_data, _maxColor, sizeof(unsigned int)*3*_width*_height);
}


PortablePixmap::~PortablePixmap() 
{
    delete[] _data;
}

void PortablePixmap::save(const std::string& filename) const 
{
    using namespace std;
    
    ofstream file;
    
    file.open(filename);

    file << "P3\n";
    file << _width << " " << _height << endl;
    file << _maxColor << endl;

    for (unsigned int i = 0; i < 3*_width*_height; i++) {
        file << _data[i] << endl;
    }

    file.close();
}

void PortablePixmap::set(unsigned int i, unsigned int j,  unsigned int r,
    unsigned int g, unsigned int b)
{
    _data[3*(i + _width*(_height - j - 1)) + 0] = r;
    _data[3*(i + _width*(_height - j - 1)) + 1] = g;
    _data[3*(i + _width*(_height - j - 1)) + 2] = b;
}

void PortablePixmap::setJET(unsigned int i, unsigned int j, float s) 
{
    float t = s*4.0f;

    float r = std::min<float>(t - 1.5f, -t + 4.5f);
    float g = std::min<float>(t - 0.5f, -t + 3.5f);
    float b = std::min<float>(t + 0.5f, -t + 2.5f);

    r = clamp(r, 0.0f, 1.0f);
    g = clamp(g, 0.0f, 1.0f);
    b = clamp(b, 0.0f, 1.0f);

    this->set(i, j, _maxColor*r, _maxColor*g, _maxColor*b);
}



float clamp(float x, float min, float max)
{
    float r = x;
    
    if (x < min) {
        r = min;
    }

    if (x > max) {
        r = max;
    }

    return r;
}