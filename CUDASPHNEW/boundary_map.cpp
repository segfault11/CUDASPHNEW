#include "util.h"
#include "boundary_map.h"
#include "portable_pixmap.h"
#include "point_in_mesh_test.h"
#include <Wm5DistPoint3Triangle3.h>
#include <cmath>
#include <iostream>

//-----------------------------------------------------------------------------
//  Aux. function declaration
//-----------------------------------------------------------------------------
void compute_distance_point_triangle(float& distance, float normal[3], 
    const float x[3], const float t0[3], const float t1[3], const float t2[3]);

//-----------------------------------------------------------------------------
//  BoundaryMap public methods
//-----------------------------------------------------------------------------
BoundaryMap::BoundaryMap(const std::string& filename)
{
    this->Load(filename);
}
//-----------------------------------------------------------------------------
BoundaryMap::BoundaryMap(const BoundaryMapConfiguration& c):
    mDx(c.dx), mCompactSupport(c.compactSupport),
        mRestDistance(c.restDistance), mState(ALLOCATED) 
{
    mMaxDist = std::max<float>(mRestDistance, mCompactSupport);    
}
//-----------------------------------------------------------------------------
BoundaryMap::~BoundaryMap() 
{
    saveDeleteArray<float>(&mNodeContentsTable);
    saveDeleteArray<unsigned int>(&mIndexMap);
    mNodeContents.GetNumCoordinates();   
    std::list<Coordinate>::const_iterator& it = mNodeContents.begin(); 
    std::list<Coordinate>::const_iterator& end = mNodeContents.end();
    Coordinate c;
    
    float* nodeContent;

    for (; it != end; it++)
    {
        mNodeContents.get(nodeContent, c);
        delete[] nodeContent;
    }
}
//-----------------------------------------------------------------------------
void BoundaryMap::AddCanvas(const TriangleMesh& mesh)
{
    if (!(mState == ALLOCATED || mState == GENERATED))
    {
        return;
    }

    //
    // compute bounding box from canvas bb
    //
    Vector3f v1 = mesh.getBoundingBox().getV1();
    Vector3f v2 = mesh.getBoundingBox().getV2();

    Vector3f diff;
    Vector3f::minus(v2, v1, diff);
    diff.scale(1.0f/mDx);
    diff.ceil();

    // Save number of grid samples in each direction
    mIMax = static_cast<unsigned int>(diff.getX()) + 1;
    mJMax = static_cast<unsigned int>(diff.getY()) + 1;
    mKMax = static_cast<unsigned int>(diff.getZ()) + 1;
    mTotalSamples = mIMax*mJMax*mKMax;

    // adjust v2 to makes sure the grid fully contains the canvas
    diff.scale(mDx);

    v2 = v1;
    v2.add(diff);

    // set bounding box
    mDomain = Rectangle3f(v1, v2);

    //
    // (re)Init signed distance field
    //
    mNodeContents.clear();
    mNodeContents.Init(mIMax, mJMax, mKMax);
    this->Dump();

    //
    //  compute signed distance field
    //
    this->computeSignedDistanceField(mNodeContents, mesh);
    mState = INITIALIZING;
}
//-----------------------------------------------------------------------------
void BoundaryMap::Generate()
{
}
//-----------------------------------------------------------------------------
void BoundaryMap::Dump() const
{
    std::cout << "Dumping Information about Boundary Map Object ... " 
        << std::endl;
    std::cout << "imax: " << mIMax << std::endl;
    std::cout << "jmax: " << mJMax << std::endl;
    std::cout << "kmax: " << mKMax << std::endl;
    std::cout << "totalSamples: " << mTotalSamples << std::endl;
    std::cout << "Size Index Map: " << mTotalSamples*4/(1024*1024) 
        << " [Mb]"<< std::endl;
    std::cout << "" << std::endl;
}
//-----------------------------------------------------------------------------
void BoundaryMap::SaveSlice(const std::string& filename) const
{
    PortablePixmap p(mIMax, mJMax, 255);
    unsigned int k = mKMax/2;
    float* nodeContent;

    /*for (unsigned int i = 0; i < mIMax; i++)
    {
        for (unsigned int j = 0; j < mJMax; j++)
        {
            if (mNodeContents.contains(Coordinate(i, j, k)))
            {
                float dist;
                mNodeContents.get(nodeContent, Coordinate(i, j, k));
                dist = nodeContent[NC_DISTANCE];

                if (dist <= 0)
                {
                    p.setJET(i, j, std::abs(dist)/mMaxDist);
                }
            }
        }
    }*/
    unsigned int idx;

    for (unsigned int i = 0; i < mIMax; i++)
    {
        for (unsigned int j = 0; j < mJMax; j++)
        {
            idx = i + mIMax*(j + mJMax*k);

            if (mIndexMap[idx] != 0)
            {
                float dist;
                dist = mNodeContentsTable[NC_NUM_ELEMENTS*mIndexMap[idx] 
                    + NC_DISTANCE];
                //dist = nodeContent[NC_DISTANCE];

                if (dist <= 0)
                {
                    p.setJET(i, j, std::abs(dist)/mMaxDist);
                }
            }
        }
    }

    p.save(filename);
}
//-----------------------------------------------------------------------------
void BoundaryMap::Save (const std::string& filename) const
{
    if (mState == ALLOCATED)
    {
        return;
    }
    
    if (mState == INITIALIZING)
    {
        std::ofstream file;
    
        file.open (filename);

        file << mRestDistance << std::endl;
        file << mCompactSupport << std::endl;
        file << mMaxDist << std::endl;
        file << mDx << std::endl;
        file << mDomain.getV1().getX() << std::endl;
        file << mDomain.getV1().getY() << std::endl;
        file << mDomain.getV1().getZ() << std::endl;
        file << mDomain.getV2().getX() << std::endl;
        file << mDomain.getV2().getY() << std::endl;
        file << mDomain.getV2().getZ() << std::endl;
        file << mIMax << std::endl;
        file << mJMax << std::endl;
        file << mKMax << std::endl;
        file << mTotalSamples << std::endl;
        file << mNodeContents.GetNumCoordinates() << std::endl;
        auto it = mNodeContents.begin();
        auto end = mNodeContents.end();

        unsigned int cntr = 0;

        for (; it != end; it++)
        {        
            float* content;
            file << (*it).i << std::endl;
            file << (*it).j << std::endl;
            file << (*it).k << std::endl;
            mNodeContents.get(content, *it); 
        
            for (unsigned int i = 0; i < NC_NUM_ELEMENTS; i++)
            {
                file << content[i] << std::endl;
            }
        
            if (cntr % 100000 == 0)
            {
                std::cout << cntr << " of " << mNodeContents.GetNumCoordinates()
                    << std::endl; 
            }
        
            cntr++;
        }
   

        file.close();
    }



}
//-----------------------------------------------------------------------------
void BoundaryMap::Load (const std::string& filename)
{
    this->Reset();
    std::ifstream file;
    file.open(filename);
    file >> mRestDistance;
    std::cout << mRestDistance << std::endl;
    file >> mCompactSupport;
    std::cout << mCompactSupport << std::endl;
    file >> mMaxDist;
    std::cout << mMaxDist << std::endl;
    file >> mDx;
    float v1[3];
    float v2[3];
    file >> v1[0];
    file >> v1[1];
    file >> v1[2];
    file >> v2[0];
    file >> v2[1];
    file >> v2[2];
    mDomain = Rectangle3f(Vector3f(v1[0], v1[1], v1[2]), 
        Vector3f(v2[0], v2[1], v2[2]));
    file >> mIMax;
    file >> mJMax;
    file >> mKMax;
    file >> mTotalSamples;
    file >> mNumCoordinates;
    mIndexMap = new unsigned int[mTotalSamples];
    mNodeContentsTable = new float[(mNumCoordinates + 1)*NC_NUM_ELEMENTS];
    mNodeContentsTable[NC_DISTANCE] = mRestDistance;
    mNodeContentsTable[NC_NORMAL_X] = 0.0f;
    mNodeContentsTable[NC_NORMAL_Y] = 0.0f;
    mNodeContentsTable[NC_NORMAL_Z] = 0.0f;


    for (unsigned int i = 0; i < mTotalSamples; i++)
    {
        mIndexMap[i] = 0;
    }

    unsigned int coord[3];
    unsigned int idx;
    float content;
    
    for (unsigned int i = 0; i < mNumCoordinates; i++)
    {
        file >> coord[0];
        file >> coord[1];
        file >> coord[2];
        idx = coord[0] + mIMax*(coord[1] + mJMax*coord[2]);
        mIndexMap[idx] = i + 1;

        for (unsigned int j = 0; j < NC_NUM_ELEMENTS; j++)
        {
            file >> content;
            mNodeContentsTable[(i + 1)*NC_NUM_ELEMENTS + j] = content;
        }
    }
    
    file.close();

    mState = GENERATED;

}
//-----------------------------------------------------------------------------
void BoundaryMap::Reset ()
{
    mNodeContents.clear();
    saveDeleteArray<unsigned int>(&mIndexMap);
    saveDeleteArray<float>(&mNodeContentsTable);
    mState = ALLOCATED;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::GetNumCoordinates () const
{
    return mNumCoordinates;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::GetNumTotalSamples () const
{
    return mTotalSamples;
}
//-----------------------------------------------------------------------------
const float* BoundaryMap::GetNodeTable () const 
{
    return mNodeContentsTable;
}
//-----------------------------------------------------------------------------
const unsigned int* BoundaryMap::GetIndexMap () const 
{
    return mIndexMap;
}
//-----------------------------------------------------------------------------
const Rectangle3f& BoundaryMap::GetDomain () const
{
    return mDomain;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::GetIMax () const
{
    return mIMax;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::GetJMax () const
{
    return mJMax;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::GetKMax() const
{
    return mKMax;
}
//-----------------------------------------------------------------------------
float BoundaryMap::GetDx() const
{
    return mDx;
}

//-----------------------------------------------------------------------------
float BoundaryMap::GetRestDistance() const
{
    return mRestDistance;
}

//-----------------------------------------------------------------------------
//  BoundaryMap private method definitions
//-----------------------------------------------------------------------------
void BoundaryMap::computeSignedDistanceField(SparseVoxelMap<float*>& map, 
        const TriangleMesh& mesh)
{
    const float* vertexList = mesh.getVertexList();
    const unsigned int* faceList = mesh.getFaceList();
    const float *v1, *v2, *v3;
    unsigned int cMin[3];
    unsigned int cMax[3];
    float x0 = mDomain.getV1().getX();
    float y0 = mDomain.getV1().getY();
    float z0 = mDomain.getV1().getZ();
    float x[3];
    float normal[3];
    float dist;

    //
    // compute distance
    //

    // for each triangle of the mesh
    for (unsigned int m = 0; m < mesh.getNumFaces(); m++) 
    {
        std::cout << m << " of " << mesh.getNumFaces() << std::endl;

        //get triangle vertices
        v1 = &vertexList[3*faceList[3*m + 0]];
        v2 = &vertexList[3*faceList[3*m + 1]];
        v3 = &vertexList[3*faceList[3*m + 2]];

        // compute bounding box of the triangle
        this->computeGridCoordinates(cMin, cMax, v1, v2, v3);

        //std::cout << "cMin[0] = " << cMin[0] << std::endl;
        //std::cout << "cMin[1] = " << cMin[1] << std::endl;
        //std::cout << "cMin[2] = " << cMin[2] << std::endl;
        //std::cout << "cMax[0] = " << cMax[0] << std::endl;
        //std::cout << "cMax[1] = " << cMax[1] << std::endl;
        //std::cout << "cMax[2] = " << cMax[2] << std::endl;

        // for each coordinate within the bounding box
        for (unsigned int k = cMin[2]; k <= cMax[2]; k++) 
        {
            for (unsigned int j = cMin[1]; j <= cMax[1]; j++) 
            {
                for (unsigned int i = cMin[0]; i <= cMax[0]; i++) 
                {
                    // translate coordinate to world coordinates
                    x[0] = x0 + i*mDx;
                    x[1] = y0 + j*mDx;
                    x[2] = z0 + k*mDx;

                    // compute distance to the triangle
                    compute_distance_point_triangle(dist, normal, x, v1, v2, v3);

                    // update sparse voxel map
                    if (dist <= mMaxDist)
                    {
                        Coordinate coord(i, j, k);

                        // if the map already contains distance information
                        // at the current coordinate
                        if (map.contains(coord)) 
                        {
                            // check if the previously computed distance is
                            // is greater than [dist].
                            float prev;
                            float* nodeContent;
                            map.get(nodeContent, coord);
                            prev = nodeContent[NC_DISTANCE];
                            if (prev > dist)
                            {
                                // if yes, update the map
                                delete[] nodeContent;
                                nodeContent = new float[NC_NUM_ELEMENTS];
                                nodeContent[NC_DISTANCE] = dist;
                                nodeContent[NC_NORMAL_X] = normal[0];
                                nodeContent[NC_NORMAL_Y] = normal[1];
                                nodeContent[NC_NORMAL_Z] = normal[2];
                                map.add(coord, nodeContent);
                            }
                        }
                        // if not yet contained in the map
                        else
                        {
                            float* nodeContent = new float[NC_NUM_ELEMENTS];
                            nodeContent[NC_DISTANCE] = dist;
                            nodeContent[NC_NORMAL_X] = normal[0];
                            nodeContent[NC_NORMAL_Y] = normal[1];
                            nodeContent[NC_NORMAL_Z] = normal[2];
                            //std::cout << dist << std::endl;
                            // just add val to the map
                            map.add(coord, nodeContent);
                        }
                    }
                }
            }
        }
    }

    //
    // compute sign of distance
    //
    std::cout << "point in mesh test" << std::endl;
    PointInMeshTest test(mesh, 20);
    
    // for each coordinate in the signed distance field
    std::list<Coordinate>::const_iterator& it = mNodeContents.begin();
    std::list<Coordinate>::const_iterator& end = mNodeContents.end();
    Coordinate c;

    for (; it != end; it++)
    { 
        // translate coordinate to world coordinate
        c = *it;
        x[0] = x0 + c.i*mDx;
        x[1] = y0 + c.j*mDx;
        x[2] = z0 + c.k*mDx;
        
        // if [x] is inside the mesh
        if (test.isContained(Vector3f(x[0], x[1], x[2])))
        {
            // negate the distances in the signed distance field
            float* nodeContent;
            mNodeContents.get(nodeContent, c);
            nodeContent[NC_DISTANCE] = -nodeContent[NC_DISTANCE]; 
        }
    }
    
    std::cout << "point in mesh test finished" << std::endl;
}
//-----------------------------------------------------------------------------
void BoundaryMap::computeGridCoordinates(unsigned int cMin[3],
        unsigned int cMax[3], const float t1[3], const float t2[3],
        const float t3[3])
{
    float min[3];
    float max[3];
    float delta = mMaxDist + mDx;

    // compute bb of triangle (t1, t2, t3)
    min[0] = std::min<float>(t1[0], std::min<float>(t2[0], t3[0]));
    min[1] = std::min<float>(t1[1], std::min<float>(t2[1], t3[1]));
    min[2] = std::min<float>(t1[2], std::min<float>(t2[2], t3[2]));
    max[0] = std::max<float>(t1[0], std::max<float>(t2[0], t3[0]));
    max[1] = std::max<float>(t1[1], std::max<float>(t2[1], t3[1]));
    max[2] = std::max<float>(t1[2], std::max<float>(t2[2], t3[2]));
    
    // extend bb to include all grid point that are closer than 
    // mMaxDist to the triangle
    min[0] -= delta;
    min[1] -= delta;
    min[2] -= delta;
    max[0] += delta;
    max[1] += delta;
    max[2] += delta;

    // clamp bb to bb of domain to not leave the domain.
    
    min[0] = std::max<float>(min[0], mDomain.getV1().getX());
    min[1] = std::max<float>(min[1], mDomain.getV1().getY());
    min[2] = std::max<float>(min[2], mDomain.getV1().getZ());
    max[0] = std::min<float>(max[0], mDomain.getV2().getX());
    max[1] = std::min<float>(max[1], mDomain.getV2().getY());
    max[2] = std::min<float>(max[2], mDomain.getV2().getZ());

    // compute min and max coordinates of the bb
    cMin[0] = static_cast<unsigned int>((min[0] - mDomain.getV1().getX())/mDx);
    cMin[1] = static_cast<unsigned int>((min[1] - mDomain.getV1().getY())/mDx);
    cMin[2] = static_cast<unsigned int>((min[2] - mDomain.getV1().getZ())/mDx);
    cMax[0] = static_cast<unsigned int>((max[0] - mDomain.getV1().getX())/mDx);
    cMax[1] = static_cast<unsigned int>((max[1] - mDomain.getV1().getY())/mDx);
    cMax[2] = static_cast<unsigned int>((max[2] - mDomain.getV1().getZ())/mDx);
}

//-----------------------------------------------------------------------------
//  Aux. function definition
//-----------------------------------------------------------------------------
void compute_distance_point_triangle(float& distance, float normal[3], 
    const float x[3], const float t0[3], const float t1[3], const float t2[3])
{
    Wm5::Vector3f point(x[0], x[1], x[2]);
    Wm5::Vector3f v0(t0[0], t0[1], t0[2]);
    Wm5::Vector3f v1(t1[0], t1[1], t1[2]);
    Wm5::Vector3f v2(t2[0], t2[1], t2[2]);
    Wm5::Triangle3f triangle(v0, v1, v2);
    Wm5::DistPoint3Triangle3f dist(point, triangle);
    distance = dist.Get();
    float a0 = dist.GetTriangleBary(0);
    float a1 = dist.GetTriangleBary(1);
    float a2 = dist.GetTriangleBary(2);
    normal[0] = x[0] - (a0*t0[0] + a1*t1[0] + a2*t2[0]);
    normal[1] = x[1] - (a0*t0[1] + a1*t1[1] + a2*t2[1]);
    normal[2] = x[2] - (a0*t0[2] + a1*t1[2] + a2*t2[2]);
    float mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + 
        normal[2]*normal[2]);
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;
}


