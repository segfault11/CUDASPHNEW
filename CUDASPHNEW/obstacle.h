#ifndef _OBSTACLE_H
#define _OBSTACLE_H

#include <string>
#include "vector3f.h"
#include "rectangle3f.h"

class ObstacleRenderer;
class ObstacleGrid;

class Obstacle
{
    friend ObstacleRenderer;
    friend ObstacleGrid;
public:
    Obstacle(const std::string& filename);
    ~Obstacle();

    void scale(const Vector3f& v);
    void scale(float s);
    void translate(const Vector3f& v);

    float computeDistance(const Vector3f& v);

    const Rectangle3f& getBoundingBox() const;
    bool isClosed() const;

private:
    inline void centralize();
private:
    float* _vertexList;
    unsigned int* _faceList;
    unsigned int _nFaces;
    unsigned int _nVertices;
    Rectangle3f _boundingBox;

};


#endif /* end of include guard: obstacle.h */