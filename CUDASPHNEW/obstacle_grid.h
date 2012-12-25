/** @file obstacle_grid.h
*** 
*** Implementation of a static boundary for sph fluid simulations as proposed 
*** in Smoothed Particle Hydrodynamics in Complex Shapes (2007)
***
*** http://individuals.iii.u-tokyo.ac.jp/~yoichiro/report/report-pdf/harada/international/2007sccg.pdf
***
**/
#ifndef _OBSTACLE_GRID_H
#define _OBSTACLE_GRID_H

#include <vector>
#include <set>
#include <list>
#include "triangle_mesh.h"

typedef struct
{
    float dx;                   /* distance between two neighbouring nodes */
    float restDistance;
    float compactSupport;
    float mass;
    float density;
}
ObstacleGridConfiguration;

class ObstacleGrid
{
    /* aux. structs */
    struct NodeTuple
    {
        float distance;         /* signed distance to boundary */
        float visWall;          /* contribution of the wall to the viscosity 
                                    force */
        float densWall;         /* contribution of the wall to the density */
        float normalDirection;  /* either 1.0f, or -1.0f */
        unsigned int normalIndex; 
    };

    struct Normal
    {
        float x, y, z;
    };

    struct GridCoordinate
    {
        unsigned int i, j, k;
        GridCoordinate(unsigned int i, unsigned int j, unsigned int k):
            i(i), j(j), k(k) {}
        ~GridCoordinate(){}
    };

public:
    ObstacleGrid(const ObstacleGridConfiguration& config);
    ~ObstacleGrid();


    void setCanvas(const TriangleMesh& obs);



    /** @brief Saves a slice of the signed distance map to a .ppm file.
    **/
    void saveDistanceMap(const std::string& filename) const;

private:
    ObstacleGrid();
    ObstacleGrid(const ObstacleGrid& orig);
    ObstacleGrid& operator = (const ObstacleGrid& orig);

    void sampleDistanceCanvas();
    void updateGridDistance(float v1[3], float v2[3], float v3[3]);


    void reset();
private:

    // general (geometric) grid information
    Rectangle3f _boundingBox;   /* Bounding box of the Obstacle Grid */ 
    unsigned int _nSamples[3];  /* number of samples in each direction */
    unsigned int _totalSamples; 
    float _dx;                  /* Grid spacing */ 
    float _compactSupport;
    float _restDistance;
    
    
    unsigned int* _nodeIndexMap; /* indices refering to node data for each 
                                  * grid point */

    const TriangleMesh* _canvas;
    std::vector<TriangleMesh*> _obstacles;
    std::set<Normal> _normal;
    std::vector<NodeTuple*> _nodeTable;
    std::list<GridCoordinate*> _coordinates;

};





#endif /* end of include guard: static_boundary.h */