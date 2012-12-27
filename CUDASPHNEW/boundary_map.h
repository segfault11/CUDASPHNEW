#ifndef _BOUNDARY_MAP_H
#define _BOUNDARY_MAP_H

#include <string>
#include "triangle_mesh.h"
#include "sparse_voxel_map.h"
#include "rectangle3f.h"

enum GenerationMode
{
    BM_NOSWAP,
    BM_HDDSWAP
};

class BoundaryMapConfiguration
{
public:
    float restDistance;
    float compactSupport;
    float dx;
};

class BoundaryMap
{
    enum
    {
        NC_DISTANCE,
        NC_NORMAL_X,
        NC_NORMAL_Y,
        NC_NORMAL_Z,
        NC_NUM_ELEMENTS
    };

    enum State
    {
        ALLOCATED,
        INITIALIZING,
        GENERATED
    };

public:
    BoundaryMap(const BoundaryMapConfiguration& c);
    ~BoundaryMap();

    void addCanvas(const TriangleMesh& mesh);
    void addObstacle(const TriangleMesh& mesh);
    void saveSlice(const std::string& filename) const;
    void generate();


    const float* getDistances() const;
    const unsigned int* getIndexMap() const;
    unsigned int getIMax() const;
    unsigned int getJMax() const;
    unsigned int getKMax() const;
    Rectangle3f getDomain() const;

    void save(const std::string& filename) const;
    void load(const std::string& filename);

    void dump() const;

protected:
    float _restDistance;
    float _compactSupport;
    float _maxDist;
    float _dx;
    Rectangle3f _domain;
    unsigned int _iMax, _jMax, _kMax;
    unsigned int _totalSamples;
    unsigned int _nCoordinates;
    
    SparseVoxelMap<float*> _nodeContents;
    unsigned int* _indexMap;
    float* _nodeContentsTable;

    State _state;

private:
    inline void computeSignedDistanceField(SparseVoxelMap<float*>& map, 
        const TriangleMesh& mesh);
    inline void computeGridCoordinates(unsigned int cMin[3],
        unsigned int cMax[3], const float t1[3], const float t2[3], 
        const float t3[3]);

};

#endif /* end of include guard: boundary_guard.h */