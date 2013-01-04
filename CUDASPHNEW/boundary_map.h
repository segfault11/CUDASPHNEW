#ifndef _BOUNDARY_MAP_H
#define _BOUNDARY_MAP_H

#include <string>
#include "triangle_mesh.h"
#include "sparse_voxel_map.h"
#include "rectangle3f.h"

enum GenerationMode
{
    GM_NOSWAP,
    GM_HDDSWAP
};

enum
{
    NC_DISTANCE,
    NC_NORMAL_X,
    NC_NORMAL_Y,
    NC_NORMAL_Z,
    NC_NUM_ELEMENTS
};

class BoundaryMapConfiguration
{
public:
    float restDistance;
    float compactSupport;
    float dx;
};
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
class BoundaryMap
{

    enum State
    {
        ALLOCATED,
        INITIALIZING,
        GENERATED
    };

    class Particle
    {
    public:
        Particle (const Vector3f& v)
            : pos(v)
        {
        }
        ~Particle ()
        {
        }
        Vector3f pos;
    };

    class ParticleManager
    {
    public:
        ParticleManager ();
        ~ParticleManager ();
        void Seed ();
    };


public:
    BoundaryMap(const std::string& filename);
    BoundaryMap(const BoundaryMapConfiguration& c);
    ~BoundaryMap();

    void AddCanvas (const TriangleMesh& mesh);
    void AddObstacle (const TriangleMesh& mesh);
    void SaveSlice (const std::string& filename) const;
    void Generate ();
    
    // Reset the objects state to ALLOCATED. Free all resources hold by the
    // object.
    void Reset ();

    // access to class
    const float* GetNodeTable () const;
    const unsigned int* GetIndexMap () const;
    unsigned int GetNumCoordinates () const;
    unsigned int GetNumTotalSamples () const;
    unsigned int GetIMax () const;
    unsigned int GetJMax () const;
    unsigned int GetKMax () const;
    float GetDx () const;
    float GetRestDistance () const;
    const Rectangle3f& GetDomain () const;

    // Save/Load boundary map
    void Save (const std::string& filename) const;
    void Load (const std::string& filename);

    void Dump () const;

protected:
    float mRestDistance;
    float mCompactSupport;
    float mMaxDist;
    float mDx;
    Rectangle3f mDomain;
    unsigned int mIMax, mJMax, mKMax;
    unsigned int mTotalSamples;
    unsigned int mNumCoordinates;
    
    SparseVoxelMap<float*> mNodeContents;
    unsigned int* mIndexMap;
    float* mNodeContentsTable;

    State mState;

private:
    inline void computeSignedDistanceField(SparseVoxelMap<float*>& map, 
        const TriangleMesh& mesh);
    inline void computeGridCoordinates(unsigned int cMin[3],
        unsigned int cMax[3], const float t1[3], const float t2[3], 
        const float t3[3]);

};

#endif /* end of include guard: boundary_guard.h */