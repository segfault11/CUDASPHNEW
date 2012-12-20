#ifndef _PORTABLE_PIXMAP_H
#define _PORTABLE_PIXMAP_H

#include <string>

class PortablePixmap
{
public:
    PortablePixmap(unsigned int width, unsigned int height, 
        unsigned int maxColor);
    ~PortablePixmap();

    void save(const std::string& filename) const;

    void set(unsigned int i, unsigned int j, 
        unsigned int r, unsigned int g, unsigned int b);

    void setJET(unsigned int i, unsigned int j, float s);

private:
    unsigned int* _data;
    unsigned int _width;
    unsigned int _height;
    unsigned int _maxColor;
};

#endif /* end of include guard: portable_pixmap.h */