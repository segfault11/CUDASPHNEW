#ifndef _BOUNDARY_MAP_H
#define _BOUNDARY_MAP_H

#include <string>
#include "triangle_mesh.h"
#include "sparse_voxel_map.h"
#include "rectangle3f.h"



class BoundaryMapConfiguration
{
public:
    float restDistance;
    float compactSupport;
    float dx;
};

class BoundaryMap
{
public:
    BoundaryMap(const BoundaryMapConfiguration& c);
    ~BoundaryMap();

    void addCanvas(const TriangleMesh& mesh);
    void addObstacle(const TriangleMesh& mesh);

    void save(const std::string& filename) const;

protected:
    SparseVoxelMap<float> _signedDistanceField;
    float _restDistance;
    float _compactSupport;
    float _maxDist;
    float _dx;
    Rectangle3f _domain;
    unsigned int _iMax, _jMax, _kMax;
    unsigned int _totalSamples;

private:
    inline void computeSignedDistanceField(SparseVoxelMap<float>& map, 
        const TriangleMesh& mesh);
    inline void computeGridCoordinates(unsigned int cMin[3],
        unsigned int cMax[3], const float t1[3], const float t2[3], 
        const float t3[3]);
};

#endif /* end of include guard: boundary_guard.h */