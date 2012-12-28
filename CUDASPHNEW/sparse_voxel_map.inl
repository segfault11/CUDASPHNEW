#include "util.h"
#include <cstdlib>
#include <limits>

#include <iostream>

//-----------------------------------------------------------------------------
template<class T>
SparseVoxelMap<T>::SparseVoxelMap(): _iMax(0), _jMax(0), _kMax(0), 
    _data(NULL), _coordinates() {}
//-----------------------------------------------------------------------------
template<class T>
SparseVoxelMap<T>::SparseVoxelMap(unsigned int iMax, unsigned int jMax, 
    unsigned int kMax): _iMax(iMax), _jMax(jMax), _kMax(kMax)
{
    _data = new T*[_iMax*_jMax*_kMax];
    memset(_data, NULL, _iMax*_jMax*_kMax*sizeof(T*));
}
//-----------------------------------------------------------------------------
template<class T>
SparseVoxelMap<T>::~SparseVoxelMap() 
{ 
    auto itBegin = _coordinates.begin();
    auto itEnd = _coordinates.end();
    Coordinate c;
    unsigned int idx;

    for (; itBegin != itEnd; itBegin++)
    {
        c = *itBegin;
        idx = c.i + _iMax*(c.j + _jMax*c.k);
        delete _data[idx];
    }

    delete[] _data;
}
//-----------------------------------------------------------------------------
template<class T>
bool SparseVoxelMap<T>::contains(const Coordinate& c) const
{
    unsigned int idx = c.i + _iMax*(c.j + _jMax*c.k);
    
    return _data[idx] != NULL;
}
//-----------------------------------------------------------------------------
template<class T>
void SparseVoxelMap<T>::add(const Coordinate& c, const T& val) 
{
    unsigned int idx = c.i + _iMax*(c.j + _jMax*c.k);

    if (this->contains(c)) 
    {
        *_data[idx] = val;
        return;
    }

    _coordinates.push_back(c);
    _data[idx] = new T(val);
}
//-----------------------------------------------------------------------------
template<class T>
bool SparseVoxelMap<T>::get(T& val, const Coordinate& c) const 
{
    if (!this->contains(c))
    {
        return false;
    }

    unsigned int idx = c.i + _iMax*(c.j + _jMax*c.k);
    
    val = *_data[idx];

    return true;
}
//-----------------------------------------------------------------------------
template<class T>
unsigned int SparseVoxelMap<T>::getNumCoordinates() const
{
    return _coordinates.size();
}
//-----------------------------------------------------------------------------
template<class T>
std::list<Coordinate>::const_iterator& SparseVoxelMap<T>::begin() const
{
    std::cout << "KARLKALSON" << std::endl;
    return _coordinates.begin();
}
//-----------------------------------------------------------------------------
template<class T>
std::list<Coordinate>::const_iterator& SparseVoxelMap<T>::end() const
{
    return _coordinates.end();
}
//-----------------------------------------------------------------------------
template<class T>
void SparseVoxelMap<T>::clear()
{
    auto itBegin = _coordinates.begin();
    auto itEnd = _coordinates.end();
    Coordinate c;
    unsigned int idx;

    for (; itBegin != itEnd; itBegin++)
    {
        c = *itBegin;
        idx = c.i + _iMax*(c.j + _jMax*c.k);
        delete _data[idx];
    }

    _coordinates.clear();

    if(_data != NULL)
    {
        delete[] _data;
    }

    _data = NULL;
}
//-----------------------------------------------------------------------------
template<class T>
void SparseVoxelMap<T>::init(unsigned int iMax, unsigned int jMax,
    unsigned int kMax)
{
    this->clear();
    _iMax = iMax;
    _jMax = jMax;
    _kMax = kMax;
    _data = new T*[_iMax*_jMax*_kMax];
    memset(_data, NULL, _iMax*_jMax*_kMax*sizeof(T*));
}
//-----------------------------------------------------------------------------
template<class T>
void SparseVoxelMap<T>::dumpCoordinates() const
{
    auto it = _coordinates.begin();
    auto end = _coordinates.end();

    for (; it != end; it++)
    {
        std::cout << (*it).i << " " << (*it).j << " " << (*it).k << std::endl;
    }


}
//-----------------------------------------------------------------------------
template<class T>
void dumpCoordinates(const SparseVoxelMap<T>& map)
{
    auto it = map.begin();
    auto end = map.end();

    for (; it != end; it++)
    {
        std::cout << (*it).i << " " << (*it).j << " " << (*it).k << std::endl;
    }

}