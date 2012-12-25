#include "boundary_map.h"
#include "portable_pixmap.h"
#include "point_in_mesh_test.h"
#include <Wm5DistPoint3Triangle3.h>
#include <cmath>
#include <iostream>


inline float compute_distance_point_triangle(const float x[3], 
    const float t0[3], const float t1[3], const float t2[3]);

BoundaryMap::BoundaryMap(const BoundaryMapConfiguration& c):
    _dx(c.dx), _compactSupport(c.compactSupport), _restDistance(c.restDistance) 
{
    _maxDist = std::max<float>(_restDistance, _compactSupport);    
}

BoundaryMap::~BoundaryMap() {}

void BoundaryMap::addCanvas(const TriangleMesh& mesh)
{
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
    _signedDistanceField.clear();
    _signedDistanceField.init(_iMax, _jMax, _kMax);

    //
    //  compute signed distance field
    //
    this->computeSignedDistanceField(_signedDistanceField, mesh);
}

void BoundaryMap::computeSignedDistanceField(SparseVoxelMap<float>& map, 
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
                    dist = compute_distance_point_triangle(x, v1, v2, v3);

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
                            
                            map.get(prev, coord);

                            if (prev > dist)
                            {
                                // if yes, update the map
                                map.add(coord, dist);
                            }
                        }
                        // if not yet contained in the map
                        else
                        {
                            //std::cout << dist << std::endl;
                            // just add val to the map
                            map.add(coord, dist);
                        }
                    }
                }
            }
        }
    }

    //
    // compute sign of distance
    //

    PointInMeshTest test(mesh, 20);
    
    // for each coordinate in the signed distance field
    std::list<Coordinate>::const_iterator& it = _signedDistanceField.begin();
    std::list<Coordinate>::const_iterator& end = _signedDistanceField.end();
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
            _signedDistanceField.get(dist, c);
            _signedDistanceField.add(c, -dist);
        }
    }
}


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

void BoundaryMap::save(const std::string& filename) const
{
    PortablePixmap p(_iMax, _jMax, 255);
    unsigned int k = _kMax/2;

    for (unsigned int i = 0; i < _iMax; i++)
    {
        for (unsigned int j = 0; j < _jMax; j++)
        {
            if (_signedDistanceField.contains(Coordinate(i, j, k)))
            {
                float dist;

                _signedDistanceField.get(dist, Coordinate(i, j, k));
                
                if (dist <= 0)
                {
                    p.setJET(i, j, std::abs(dist)/_maxDist);
                }
            }
        }
    }

    p.save(filename);
}


float compute_distance_point_triangle(const float x[3], 
    const float t0[3], const float t1[3], const float t2[3])
{
    Wm5::Vector3f point(x[0], x[1], x[2]);
    Wm5::Vector3f v0(t0[0], t0[1], t0[2]);
    Wm5::Vector3f v1(t1[0], t1[1], t1[2]);
    Wm5::Vector3f v2(t2[0], t2[1], t2[2]);
    Wm5::Triangle3f triangle(v0, v1, v2);

    return Wm5::DistPoint3Triangle3f(point, triangle).Get();
}
