#ifndef _POINT_IN_MESH_TEST_H
#define _POINT_IN_MESH_TEST_H

#include "triangle_mesh.h"
#include "rectangle3f.h"
#include "vector3f.h"
#include <Wm5IntrRay3Triangle3.h>
#include <list>
#include <string>

class PointInMeshTest
{
public:
    PointInMeshTest(const TriangleMesh& mesh);
    PointInMeshTest(const TriangleMesh& mesh, unsigned int nSamples);
    ~PointInMeshTest();

    bool isContained(const Vector3f& point) const;
    bool isInside(const Vector3f& point) const;

    void save(const std::string& filename) const;

private:
    inline void hashTriangles();
    inline bool test(bool& isInside, const Wm5::Ray3f& ray, 
        const std::list<unsigned int>& triangeList) const;

private:
    const TriangleMesh* _mesh;
    Rectangle3f _boundingBox;
    float _dx;
    float _dy;
    unsigned int _nSamples;
    std::list<unsigned int>* _hashGrid;
    unsigned int _nMaxTrials;
};

#endif /* end of include guard: point_in_mesh_test.h */