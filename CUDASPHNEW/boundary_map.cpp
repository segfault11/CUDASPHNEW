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
BoundaryMap::BoundaryMap()
{
 //   memset(this, 0, sizeof(BoundaryMap));
}
//-----------------------------------------------------------------------------
BoundaryMap::BoundaryMap(const BoundaryMapConfiguration& c):
    _dx(c.dx), _compactSupport(c.compactSupport),
        _restDistance(c.restDistance), _state(ALLOCATED) 
{
    _maxDist = std::max<float>(_restDistance, _compactSupport);    
}
//-----------------------------------------------------------------------------
BoundaryMap::~BoundaryMap() 
{
    saveDeleteArray<float>(&_nodeContentsTable);
    saveDeleteArray<unsigned int>(&_indexMap);
    _nodeContents.getNumCoordinates();   
    std::list<Coordinate>::const_iterator& it = _nodeContents.begin(); 
    std::list<Coordinate>::const_iterator& end = _nodeContents.end();
    Coordinate c;
    
    float* nodeContent;

    for (; it != end; it++)
    {
        _nodeContents.get(nodeContent, c);
        delete[] nodeContent;
    }
}
//-----------------------------------------------------------------------------
void BoundaryMap::addCanvas(const TriangleMesh& mesh)
{
    if (!(_state == ALLOCATED || _state == GENERATED))
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
    diff.scale(1.0f/_dx);
    diff.ceil();

    // save number of grid samples in each direction
    _iMax = static_cast<unsigned int>(diff.getX()) + 1;
    _jMax = static_cast<unsigned int>(diff.getY()) + 1;
    _kMax = static_cast<unsigned int>(diff.getZ()) + 1;
    _totalSamples = _iMax*_jMax*_kMax;

    // adjust v2 to makes sure the grid fully contains the canvas
    diff.scale(_dx);

    v2 = v1;
    v2.add(diff);

    // set bounding box
    _domain = Rectangle3f(v1, v2);

    //
    // (re)init signed distance field
    //
    _nodeContents.clear();
    _nodeContents.init(_iMax, _jMax, _kMax);
    this->dump();

    //
    //  compute signed distance field
    //
    this->computeSignedDistanceField(_nodeContents, mesh);

    _state = INITIALIZING;
}
//-----------------------------------------------------------------------------
void BoundaryMap::generate()
{
}
//-----------------------------------------------------------------------------
void BoundaryMap::dump() const
{
    std::cout << "Dumping Information about Boundary Map Object ... " 
        << std::endl;
    std::cout << "imax: " << _iMax << std::endl;
    std::cout << "jmax: " << _jMax << std::endl;
    std::cout << "kmax: " << _kMax << std::endl;
    std::cout << "totalSamples: " << _totalSamples << std::endl;
    std::cout << "Size Index Map: " << _totalSamples*4/(1024*1024) 
        << " [Mb]"<< std::endl;
    std::cout << "" << std::endl;
}
//-----------------------------------------------------------------------------
void BoundaryMap::saveSlice(const std::string& filename) const
{
    PortablePixmap p(_iMax, _jMax, 255);
    unsigned int k = _kMax/2;
    float* nodeContent;

    /*for (unsigned int i = 0; i < _iMax; i++)
    {
        for (unsigned int j = 0; j < _jMax; j++)
        {
            if (_nodeContents.contains(Coordinate(i, j, k)))
            {
                float dist;
                _nodeContents.get(nodeContent, Coordinate(i, j, k));
                dist = nodeContent[NC_DISTANCE];

                if (dist <= 0)
                {
                    p.setJET(i, j, std::abs(dist)/_maxDist);
                }
            }
        }
    }*/
    unsigned int idx;

    for (unsigned int i = 0; i < _iMax; i++)
    {
        for (unsigned int j = 0; j < _jMax; j++)
        {
            idx = i + _iMax*(j + _jMax*k);

            if (_indexMap[idx] != 0)
            {
                float dist;
                dist = _nodeContentsTable[NC_NUM_ELEMENTS*_indexMap[idx] 
                    + NC_DISTANCE];
                //dist = nodeContent[NC_DISTANCE];

                if (dist <= 0)
                {
                    p.setJET(i, j, std::abs(dist)/_maxDist);
                }
            }
        }
    }

    p.save(filename);
}
//-----------------------------------------------------------------------------
void BoundaryMap::save(const std::string& filename) const
{
    std::ofstream file;
    
    file.open(filename);

    file << _restDistance << std::endl;
    file << _compactSupport << std::endl;
    file << _maxDist << std::endl;
    file << _dx << std::endl;
    file << _domain.getV1().getX() << std::endl;
    file << _domain.getV1().getY() << std::endl;
    file << _domain.getV1().getZ() << std::endl;
    file << _domain.getV2().getX() << std::endl;
    file << _domain.getV2().getY() << std::endl;
    file << _domain.getV2().getZ() << std::endl;
    file << _iMax << std::endl;
    file << _jMax << std::endl;
    file << _kMax << std::endl;
    file << _totalSamples << std::endl;
    file << _nodeContents.getNumCoordinates() << std::endl;
    
    auto it = _nodeContents.begin();
    auto end = _nodeContents.end();

    unsigned int cntr = 0;

    for (; it != end; it++)
    {        
        float* content;
        file << (*it).i << std::endl;
        file << (*it).j << std::endl;
        file << (*it).k << std::endl;
        _nodeContents.get(content, *it); 
        
        for (unsigned int i = 0; i < NC_NUM_ELEMENTS; i++)
        {
            file << content[i] << std::endl;
        }
        
        if (cntr % 100000 == 0)
        {
            std::cout << cntr << " of " << _nodeContents.getNumCoordinates() << std::endl; 
        }
        
        cntr++;
    }
   

    file.close();

}
//-----------------------------------------------------------------------------
void BoundaryMap::load(const std::string& filename)
{
    std::ifstream file;
    file.open(filename);
    file >> _restDistance;
    std::cout << _restDistance << std::endl;
    file >> _compactSupport;
    std::cout << _compactSupport << std::endl;
    file >> _maxDist;
    std::cout << _maxDist << std::endl;
    file >> _dx;
    float v1[3];
    float v2[3];
    file >> v1[0];
    file >> v1[1];
    file >> v1[2];
    file >> v2[0];
    file >> v2[1];
    file >> v2[2];
    _domain = Rectangle3f(Vector3f(v1[0], v1[1], v1[2]), 
        Vector3f(v2[0], v2[1], v2[2]));
    file >> _iMax;
    file >> _jMax;
    file >> _kMax;
    file >> _totalSamples;
    file >> _nCoordinates;
    _indexMap = new unsigned int[_totalSamples];
    _nodeContentsTable = new float[(_nCoordinates + 1)*NC_NUM_ELEMENTS];
    _nodeContentsTable[NC_DISTANCE] = _restDistance;

    for (unsigned int i = 0; i < _totalSamples; i++)
    {
        _indexMap[i] = 0;
    }

    unsigned int coord[3];
    unsigned int idx;
    float content;
    
    for (unsigned int i = 0; i < _nCoordinates; i++)
    {
        file >> coord[0];
        file >> coord[1];
        file >> coord[2];
        idx = coord[0] + _iMax*(coord[1] + _jMax*coord[2]);
        _indexMap[idx] = i + 1;

        for (unsigned int j = 0; j < NC_NUM_ELEMENTS; j++)
        {
            file >> content;
            _nodeContentsTable[(i + 1)*NC_NUM_ELEMENTS + j] = content;
        }
    }
    

    file.close();

}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::getNumCoordinates() const
{
    return _nCoordinates;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::getNumTotalSamples() const
{
    return _totalSamples;
}
//-----------------------------------------------------------------------------
const float* BoundaryMap::getNodeTable() const 
{
    return _nodeContentsTable;
}
//-----------------------------------------------------------------------------
const unsigned int* BoundaryMap::getIndexMap() const 
{
    return _indexMap;
}
//-----------------------------------------------------------------------------
const Rectangle3f& BoundaryMap::getDomain() const
{
    return _domain;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::getIMax() const
{
    return _iMax;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::getJMax() const
{
    return _jMax;
}
//-----------------------------------------------------------------------------
unsigned int BoundaryMap::getKMax() const
{
    return _kMax;
}
//-----------------------------------------------------------------------------
float BoundaryMap::getDx() const
{
    return _dx;
}

//-----------------------------------------------------------------------------
float BoundaryMap::getRestDistance() const
{
    return _restDistance;
}

//-----------------------------------------------------------------------------
//  BoundaryMap private methods
//-----------------------------------------------------------------------------
void BoundaryMap::computeSignedDistanceField(SparseVoxelMap<float*>& map, 
        const TriangleMesh& mesh)
{
    const float* vertexList = mesh.getVertexList();
    const unsigned int* faceList = mesh.getFaceList();
    const float *v1, *v2, *v3;
    unsigned int cMin[3];
    unsigned int cMax[3];
    float x0 = _domain.getV1().getX();
    float y0 = _domain.getV1().getY();
    float z0 = _domain.getV1().getZ();
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
                    x[0] = x0 + i*_dx;
                    x[1] = y0 + j*_dx;
                    x[2] = z0 + k*_dx;

                    // compute distance to the triangle
                    compute_distance_point_triangle(dist, normal, x, v1, v2, v3);

                    // update sparse voxel map
                    if (dist <= _maxDist)
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
    std::list<Coordinate>::const_iterator& it = _nodeContents.begin();
    std::list<Coordinate>::const_iterator& end = _nodeContents.end();
    Coordinate c;

    for (; it != end; it++)
    { 
        // translate coordinate to world coordinate
        c = *it;
        x[0] = x0 + c.i*_dx;
        x[1] = y0 + c.j*_dx;
        x[2] = z0 + c.k*_dx;
        
        // if [x] is inside the mesh
        if (test.isContained(Vector3f(x[0], x[1], x[2])))
        {
            // negate the distances in the signed distance field
            float* nodeContent;
            _nodeContents.get(nodeContent, c);
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
    float delta = _maxDist + _dx;

    // compute bb of triangle (t1, t2, t3)
    min[0] = std::min<float>(t1[0], std::min<float>(t2[0], t3[0]));
    min[1] = std::min<float>(t1[1], std::min<float>(t2[1], t3[1]));
    min[2] = std::min<float>(t1[2], std::min<float>(t2[2], t3[2]));
    max[0] = std::max<float>(t1[0], std::max<float>(t2[0], t3[0]));
    max[1] = std::max<float>(t1[1], std::max<float>(t2[1], t3[1]));
    max[2] = std::max<float>(t1[2], std::max<float>(t2[2], t3[2]));
    
    // extend bb to include all grid point that are closer than 
    // _maxDist to the triangle
    min[0] -= delta;
    min[1] -= delta;
    min[2] -= delta;
    max[0] += delta;
    max[1] += delta;
    max[2] += delta;

    // clamp bb to bb of domain to not leave the domain.
    
    min[0] = std::max<float>(min[0], _domain.getV1().getX());
    min[1] = std::max<float>(min[1], _domain.getV1().getY());
    min[2] = std::max<float>(min[2], _domain.getV1().getZ());
    max[0] = std::min<float>(max[0], _domain.getV2().getX());
    max[1] = std::min<float>(max[1], _domain.getV2().getY());
    max[2] = std::min<float>(max[2], _domain.getV2().getZ());

    // compute min and max coordinates of the bb
    cMin[0] = static_cast<unsigned int>((min[0] - _domain.getV1().getX())/_dx);
    cMin[1] = static_cast<unsigned int>((min[1] - _domain.getV1().getY())/_dx);
    cMin[2] = static_cast<unsigned int>((min[2] - _domain.getV1().getZ())/_dx);
    cMax[0] = static_cast<unsigned int>((max[0] - _domain.getV1().getX())/_dx);
    cMax[1] = static_cast<unsigned int>((max[1] - _domain.getV1().getY())/_dx);
    cMax[2] = static_cast<unsigned int>((max[2] - _domain.getV1().getZ())/_dx);
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


