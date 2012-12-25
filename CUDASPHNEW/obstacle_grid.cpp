#include "obstacle_grid.h"
#include <Wm5DistPoint3Triangle3.h>
#include "util.h"
#include <limits>
#include "portable_pixmap.h"
#include "point_in_mesh_test.h"
#include <cmath>
#include "cgtk\include\clock.h"

inline float minf(float f0, float f1);
inline float maxf(float f0, float f1);
inline int mini(int f0, int f1);
inline int maxi(int f0, int f1);


inline float compute_distance_point_triangle(float x[3], float t0[3], 
    float t1[3], float t2[3]);

ObstacleGrid::ObstacleGrid(const ObstacleGridConfiguration& config):
    _dx(config.dx), _compactSupport(config.compactSupport),
    _restDistance(config.restDistance)
{
    // first node tuple in the node table refers to node with infinity 
    // distance to the obstacles/canvas
    #undef max;

    float inf = std::numeric_limits<float>::max();

    NodeTuple* n = new NodeTuple();
    n->distance = inf;

    _nodeTable.push_back(n);

}

ObstacleGrid::~ObstacleGrid()
{

}


void ObstacleGrid::setCanvas(const TriangleMesh& obs)
{
    // TODO:
    // reset to initial state ...
    // check if canvas contains all obstacles (if obstacles were already added)

    // set canvas
    _canvas = &obs;

    // compute bounding box from canvas bb
    Vector3f v1 = _canvas->getBoundingBox().getV1();
    Vector3f v2 = _canvas->getBoundingBox().getV2();

    Vector3f diff;
    Vector3f::minus(v2, v1, diff);
    diff.scale(1.0f/_dx);
    diff.ceil();

    // save number of grid samples in each direction
    _nSamples[0] = static_cast<unsigned int>(diff.getX());
    _nSamples[1] = static_cast<unsigned int>(diff.getY());
    _nSamples[2] = static_cast<unsigned int>(diff.getZ());
    _totalSamples = _nSamples[0]*_nSamples[1]*_nSamples[2];

    // adjust v2 to makes sure the grid fully contains the canvas
    diff.scale(_dx);

    v2 = v1;
    v2.add(diff);

    // set bounding box
    _boundingBox = Rectangle3f(v1, v2);

    obs.getBoundingBox().dump();
    _boundingBox.dump();
    std::cout << _totalSamples << std::endl;

    // allocate memory for the node index map
    _nodeIndexMap = new unsigned int[_totalSamples];
    memset(_nodeIndexMap, 0, _totalSamples*sizeof(unsigned int));

    // sample signed distance field to the canvas 
    this->sampleDistanceCanvas();


}

void ObstacleGrid::sampleDistanceCanvas()
{
    // for each triangle of the canvas ...
    for (unsigned int i = 0; i < _canvas->_nFaces; i++) 
    {
        // update regions in the grid which lie within the bounding box of the triangle
        updateGridDistance(&_canvas->_vertexList[3*_canvas->_faceList[3*i + 0]],
            &_canvas->_vertexList[3*_canvas->_faceList[3*i + 1]],
            &_canvas->_vertexList[3*_canvas->_faceList[3*i + 2]]);
        std::cout << i << " of " << _canvas->_nFaces << " faces done" << std::endl;
    }


    PointInMeshTest test(*_canvas, 10);
    test.save("pittest.ppm");

    // compute sign of each distance
    float gridWC[3];  // world coord. for each grid index
    int idx; // index for node index map
    NodeTuple* t;
    

    // for all coordinates within a narrow band of the surface
    // compute the sign of the distance
    std::list<GridCoordinate*>::iterator it = _coordinates.begin();
    std::list<GridCoordinate*>::iterator end = _coordinates.end();
    unsigned int cntr = 0;

    for ( ; it != end; it++) 
    {
        if (cntr % 100000 == 0)
        {
            std::cout << cntr << "/" << _coordinates.size() << std::endl;    
        }

        cntr++;

        gridWC[0] = _boundingBox.getV1().getX() + _dx*(*it)->i; 
        gridWC[1] = _boundingBox.getV1().getY() + _dx*(*it)->j; 
        gridWC[2] = _boundingBox.getV1().getZ() + _dx*(*it)->k; 
        idx = (*it)->i + _nSamples[0]*((*it)->j + (*it)->k*_nSamples[1]);
        

        //cgtkClockStart();
        if (test.isContained(Vector3f(gridWC[0], gridWC[1], gridWC[2])))
        {
            t = _nodeTable[_nodeIndexMap[idx]];
            t->distance = -t->distance;
        }
        //cgtkClockDumpElapsed();
    }


}

void ObstacleGrid::updateGridDistance(float v1[3], float v2[3], float v3[3])
{
    // compute bounding box of the triangle
    float min[3];
    float max[3];

    min[0] = minf(minf(v1[0], v2[0]), v3[0]);
    min[1] = minf(minf(v1[1], v2[1]), v3[1]);
    min[2] = minf(minf(v1[2], v2[2]), v3[2]);

    max[0] = maxf(maxf(v1[0], v2[0]), v3[0]);
    max[1] = maxf(maxf(v1[1], v2[1]), v3[1]);
    max[2] = maxf(maxf(v1[2], v2[2]), v3[2]);

    // expand bb with respect to h and dx
    float d = maxf(_compactSupport, _restDistance) + _dx;
    
    min[0] -= d;
    min[1] -= d;
    min[2] -= d;

    max[0] += d;
    max[1] += d;
    max[2] += d;

    // compute grid indices for bounding box 
    // NOTE: the bounding box contains all grid indices
    int idxMin[3];
    int idxMax[3];

    idxMin[0] = static_cast<int>((min[0] - 
        _boundingBox.getV1().getX())/_dx) + 1;
    idxMin[1] = static_cast<int>((min[1] - 
        _boundingBox.getV1().getY())/_dx) + 1;
    idxMin[2] = static_cast<int>((min[2] - 
        _boundingBox.getV1().getZ())/_dx) + 1;

    idxMax[0] = static_cast<int>((max[0] - 
        _boundingBox.getV1().getX())/_dx);
    idxMax[1] = static_cast<int>((max[1] - 
        _boundingBox.getV1().getY())/_dx);
    idxMax[2] = static_cast<int>((max[2] - 
        _boundingBox.getV1().getZ())/_dx);

    // clamp indices, if they lie outside
    idxMin[0] = maxi(idxMin[0], 0);
    idxMin[1] = maxi(idxMin[1], 0);
    idxMin[2] = maxi(idxMin[2], 0);
    idxMax[0] = mini(idxMax[0], _nSamples[0] - 1);
    idxMax[1] = mini(idxMax[1], _nSamples[1] - 1);
    idxMax[2] = mini(idxMax[2], _nSamples[2] - 1);

    // for each grid point compute its distance to the triangle
    float gridWC[3];  // world coord. for each grid index
    float dist;      
    int idx; // index for node index map
    NodeTuple* t;

    // compute distance
    for (int k = idxMin[2]; k <= idxMax[2]; k++) {
        for (int j = idxMin[1]; j <= idxMax[1]; j++) {
            for (int i = idxMin[0]; i <= idxMax[0]; i++) {

                // translate grid indices to world coordinates
                gridWC[0] = _boundingBox.getV1().getX() + _dx*i; 
                gridWC[1] = _boundingBox.getV1().getY() + _dx*j; 
                gridWC[2] = _boundingBox.getV1().getZ() + _dx*k; 
  
                // compute distance to triangle
                dist = compute_distance_point_triangle(gridWC, v1, v2, v3);

                idx = i + _nSamples[0]*(j + k*_nSamples[1]);
            
                // only grid nodes that are closer to the obstacle than the rest distance or 
                // the compact support are interesting...
                if (dist <= d) {

                    // if no "real" entry for this grid node exists ...
                    if (_nodeIndexMap[idx] == 0) {

                        // add a new entry to the node table and update the node index map
                        t = new NodeTuple();

                        t->distance = dist;
                        _nodeTable.push_back(t);
                        _nodeIndexMap[idx] = _nodeTable.size() - 1;

                        // save not empty coordinates
                        _coordinates.push_back(new GridCoordinate(i,j,k)); 

                    }  else {

                        // else update table if dist is smaller
                        t = _nodeTable[_nodeIndexMap[idx]];

                        if (t->distance > dist) {
                            t->distance = dist;
                        }
                    }
                }
            }
        }
    }

}

void ObstacleGrid::saveDistanceMap(const std::string& filename) const
{
    unsigned int k = _nSamples[2]/2;
    unsigned int idx;
    float maxDist = _dx + std::max<float>(_compactSupport, _restDistance);
    float d;
    PortablePixmap p(_nSamples[0], _nSamples[1], 255);

    for (unsigned int i = 0; i < _nSamples[0]; i++) {
        for (unsigned int j = 0; j < _nSamples[1]; j++) {
            idx = i + _nSamples[0]*(j + k*_nSamples[1]);
            
            d = _nodeTable[_nodeIndexMap[idx]]->distance;

            if (std::abs(d) <= maxDist) {
                if (d <= 0.0) {
                    p.set(i, j, 255, 0, 0);
                } else {
                    p.set(i, j, 0, 255, 0);
                }
                
                
                //p.setJET(i, j, d/maxDist);
            }

        }    
    }

    p.save(filename);
}

float minf(float f0, float f1) {
    return f0 < f1 ? f0 : f1;
}

float maxf(float f0, float f1) {
    return f0 > f1 ? f0 : f1;
}

int mini(int f0, int f1) {
    return f0 < f1 ? f0 : f1;
}

int maxi(int f0, int f1) {
    return f0 > f1 ? f0 : f1;
}

float compute_distance_point_triangle(float x[3], float t0[3], 
    float t1[3], float t2[3])
{
    Wm5::Vector3f point(x[0], x[1], x[2]);
    Wm5::Vector3f v0(t0[0], t0[1], t0[2]);
    Wm5::Vector3f v1(t1[0], t1[1], t1[2]);
    Wm5::Vector3f v2(t2[0], t2[1], t2[2]);
    Wm5::Triangle3f triangle(v0, v1, v2);
    return Wm5::DistPoint3Triangle3f(point, triangle).Get();
}
