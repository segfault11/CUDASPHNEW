#include "point_in_mesh_test.h"
#include "vector3f.h"
#include "portable_pixmap.h"
#include <iostream>
#include <cmath>

#include "cgtk\include\clock.h"

PointInMeshTest::PointInMeshTest(const TriangleMesh& mesh, 
    unsigned int nSamples): _mesh(&mesh), _nMaxTrials(5), _nSamples(nSamples), 
        _boundingBox(mesh.getBoundingBox())
{
    Vector3f diff;

    Vector3f::minus(_boundingBox.getV2(), _boundingBox.getV1(), diff);

    _dx = diff.getX()/(_nSamples - 1);
    _dy = diff.getY()/(_nSamples - 1);

    _hashGrid = new std::list<unsigned int>[(_nSamples)*(_nSamples)];

    this->hashTriangles();
}

PointInMeshTest::~PointInMeshTest() 
{
    delete[] _hashGrid;
}


bool PointInMeshTest::isContained(const Vector3f& point) const
{
    if (!_boundingBox.contains(point)) {
        return false;
    }

    std::srand(std::time(NULL));
    

    float xs = _boundingBox.getV1().getX();
    float ys = _boundingBox.getV1().getY();

    unsigned int i, j, idx;
    
    i = static_cast<unsigned int>((point.getX() - xs)/_dx);
    j = static_cast<unsigned int>((point.getY() - ys)/_dy);
    idx = i + _nSamples*j;

    if (idx >= _nSamples*_nSamples) 
    {
        std::cout << "Warning idx is " << idx  << std::endl;
        system("pause");
    }


    std::list<unsigned int>& triangleList = _hashGrid[idx];

    unsigned int face;
    float  *v1, *v2, *v3;
    unsigned int cnt = 0;




    bool succ = true;

    do {
        std::list<unsigned int>::const_iterator& it = triangleList.begin();
        std::list<unsigned int>::const_iterator& end = triangleList.end();
        cnt = 0;

        float off= static_cast<float>(std::rand() % 999)/1000.0f;

        float rayX = off*_dx;
        float rayY = off*_dy;
        float rayZ = _boundingBox.getV2().getZ() - _boundingBox.getV1().getZ();

        Wm5::Vector3f origin(point.getX(), point.getY(), point.getZ());
        Wm5::Ray3f ray(origin, Wm5::Vector3f(rayX, rayX, rayZ));



        for ( ; it != end; it++) 
        {
            face = *it;
            v1 = &_mesh->_vertexList[3*_mesh->_faceList[3*face + 0]];
            v2 = &_mesh->_vertexList[3*_mesh->_faceList[3*face + 1]];
            v3 = &_mesh->_vertexList[3*_mesh->_faceList[3*face + 2]];
    
            Wm5::Triangle3f tr(Wm5::Vector3f(v1[0], v1[1], v1[2]),
                Wm5::Vector3f(v2[0], v2[1], v2[2]),
                Wm5::Vector3f(v3[0], v3[1], v3[2]));

            Wm5::IntrRay3Triangle3f testtt(ray, tr);

            if (testtt.Find())
            {

                if (testtt.GetTriBary0() == 0.0f || testtt.GetTriBary1() == 0.0f
                    || testtt.GetTriBary2() == 0.0f) 
                {
                    std::cout << "Dies ist eine Warnung" << std::endl;
                    succ = false;
                }
                else
                {
                    succ = true;
                }

           
                cnt++;
            }
        
        }
    }
    while (!succ);

  
    if (cnt % 2 == 1) 
    {
        return true;
    } 
    else 
    {
        return false;
    }
}

bool PointInMeshTest::isInside(const Vector3f& point) const
{
    float  *v1, *v2, *v3;
    unsigned int cnt = 0;

    Wm5::Vector3f origin(point.getX(), point.getY(), point.getZ());
    Wm5::Ray3f ray(origin, Wm5::Vector3f(1.0f, 0.0f, 0.0f));

 
    for (unsigned int i = 0; i < _mesh->_nFaces; i++)
    {
        v1 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 0]];
        v2 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 1]];
        v3 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 2]];
        
        Wm5::Triangle3f tr(Wm5::Vector3f(v1[0], v1[1], v1[2]),
            Wm5::Vector3f(v2[0], v2[1], v2[2]),
            Wm5::Vector3f(v3[0], v3[1], v3[2]));


        Wm5::IntrRay3Triangle3f test(ray, tr);

        if (test.Test())
        {
            cnt++;
        }
    }

    if (cnt % 2 == 1) 
    {
        return true;
    } 
    else 
    {
        return false;
    }

}

void PointInMeshTest::hashTriangles()
{
    float  *v1, *v2, *v3;
    float bbMin[2];
    float bbMax[2];
    float x, y;
    float xs = _boundingBox.getV1().getX();
    float ys = _boundingBox.getV1().getY();
    unsigned int ii, jj;
    unsigned int iiMax, jjMax;
    unsigned int idx;

    for (unsigned int i = 0; i < _mesh->_nFaces; i++) 
    {
        v1 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 0]];
        v2 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 1]];
        v3 = &_mesh->_vertexList[3*_mesh->_faceList[3*i + 2]];

        bbMin[0] = std::min<float>(v1[0], std::min<float>(v2[0], v3[0]));
        bbMin[1] = std::min<float>(v1[1], std::min<float>(v2[1], v3[1]));
        bbMax[0] = std::max<float>(v1[0], std::max<float>(v2[0], v3[0]));
        bbMax[1] = std::max<float>(v1[1], std::max<float>(v2[1], v3[1]));

        ii = static_cast<unsigned int>((bbMin[0] - xs)/_dx);
        jj = static_cast<unsigned int>((bbMin[1] - ys)/_dy);
        iiMax = static_cast<unsigned int>((bbMax[0] - xs)/_dx);
        jjMax = static_cast<unsigned int>((bbMax[1] - ys)/_dy);

        for (unsigned int u = ii; u <= iiMax; u++) 
        {
            for (unsigned int v = jj; v <= jjMax; v++) 
            {
                idx = u + _nSamples*v;
                _hashGrid[idx].push_back(i);
            }    
        }

    }
}

void PointInMeshTest::save(const std::string& filename) const
{
    PortablePixmap p(_nSamples, _nSamples, 255);

    for (unsigned int i = 0; i < _nSamples; i++) {
        for (unsigned int j = 0; j < _nSamples; j++) {
            if (0 != _hashGrid[i + _nSamples*j].size()) {
                p.set(i, j, 255, 0, 0);
            }
        }
    }

    p.save(filename);
}