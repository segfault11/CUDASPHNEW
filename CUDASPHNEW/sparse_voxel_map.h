#ifndef _SPARSE_VOXEL_MAP_H
#define _SPARSE_VOXEL_MAP_H

#include <list>
#include <vector>

/** @class Coordinate
**/
class Coordinate
{
public:
    unsigned int i, j, k;

public:
    Coordinate(): i(0), j(0), k(0) {}

    Coordinate(unsigned int i, unsigned int j, unsigned int k):
        i(i), j(j), k(k) {}

    Coordinate& operator = (const Coordinate& orig) 
    {
        i = orig.i;    
        j = orig.j;    
        k = orig.k;    
    
        return *this;
    }
    
    ~Coordinate() {}
};

/** @class SparseVoxelMap
**/
template<class T>
class SparseVoxelMap
{
public:
    SparseVoxelMap();
    SparseVoxelMap(unsigned int iMax, unsigned int jMax, unsigned int kMax);
    ~SparseVoxelMap();

    /** @brief Adds a datum [val] at coordinate position [c] to the sparse \
    ***        voxel map.
    ***
    *** Note: If a datum has prevously been inserted at [c] it will be 
    *** overwritten with the new datum.
    ***
    *** @param c Position in the map where the datum should be added.
    *** @param val Datum to be added.
    **/
    void add(const Coordinate& c, const T& val);

    /** @brief Checks if a datum is contained at coordinate [c].
    **/
    bool contains(const Coordinate& c) const;
    bool get(T& val, const Coordinate& c) const;


    void init(unsigned int iMax, unsigned int jMax, unsigned int kMax);
    void clear();


    std::list<Coordinate>::const_iterator& begin() const;
    std::list<Coordinate>::const_iterator& end() const;

private:
    unsigned int _iMax;
    unsigned int _jMax;
    unsigned int _kMax;
    std::list<Coordinate> _coordinates;
    T** _data;
};

#include "sparse_voxel_map.inl"

#endif /* end of include guard: sparse_voxel_map.h */