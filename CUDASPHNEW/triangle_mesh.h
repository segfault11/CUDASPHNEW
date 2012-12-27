#ifndef _OBSTACLE_H
#define _OBSTACLE_H

#include <string>
#include "vector3f.h"
#include "rectangle3f.h"

class ObstacleRenderer;
class ObstacleGrid;
class PointInMeshTest;

class TriangleMesh
{
    friend PointInMeshTest;
    friend ObstacleRenderer;
    friend ObstacleGrid;
public:
    TriangleMesh(const std::string& filename);
    ~TriangleMesh();

    void scale(const Vector3f& v);
    void scale(float s);
    void translate(const Vector3f& v);
    void fit(const Rectangle3f& rect);

    const Rectangle3f& getBoundingBox() const;
    bool isClosed() const;

    unsigned int getNumFaces() const;
    unsigned int getNumVertices() const;
    const unsigned int* getFaceList() const;
    const float* getVertexList() const;

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