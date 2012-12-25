#include "triangle_mesh.h"
#include "cgtk\include\geometry.h"
#include "util.h"
#include <limits>
#include <Wm5DistPoint3Triangle3.h>

inline Rectangle3f compute_bounding_box(float* vertexList, 
    unsigned int nVertices);

TriangleMesh::TriangleMesh(const std::string& filename)
{
    CGTKObjFile* obs = cgtkObjFileAlloc(filename.c_str());
    
    if (obs == NULL) {
        THROW_EXCEPTION("Could not load obstacle from file.");
    }

    _nFaces = obs->nFaces;
    _nVertices = obs->nVertices;

    _vertexList = new float[_nVertices*3];
    _faceList = new unsigned int[_nFaces*3];

    memcpy(_vertexList, obs->vertices, sizeof(float)*3*obs->nVertices);
    memcpy(_faceList, obs->indices, sizeof(unsigned int)*3*obs->nFaces);

    _boundingBox = compute_bounding_box(_vertexList, _nVertices);

    this->centralize();


    cgtkObjFileFree(&obs);
}

TriangleMesh::~TriangleMesh()
{
    delete[] _vertexList;
    delete[] _faceList;
}

const Rectangle3f& TriangleMesh::getBoundingBox() const
{
    return _boundingBox;
}

void TriangleMesh::translate(const Vector3f& v)
{
    for (unsigned int i = 0; i < _nVertices; i++) {
        _vertexList[3*i + 0] += v.getX();
        _vertexList[3*i + 1] += v.getY();
        _vertexList[3*i + 2] += v.getZ();
    }
}

void TriangleMesh::scale(float s)
{
    // translate object around origin
    Vector3f mid;
    Vector3f::minus(_boundingBox.getV2(), _boundingBox.getV1(), mid);
    mid.scale(0.5f);
    Vector3f::plus(_boundingBox.getV1(), mid, mid);
    mid.scale(-1.0f);

    this->translate(mid);

    _boundingBox = compute_bounding_box(_vertexList, _nVertices);

    // scale object
    _boundingBox.scale(s);
    for (unsigned int i = 0; i < _nVertices; i++) {
        _vertexList[3*i + 0] *= s;
        _vertexList[3*i + 1] *= s;
        _vertexList[3*i + 2] *= s;
    }



    // translate object back
    mid.scale(-1.0f);

    this->translate(mid);

    _boundingBox = compute_bounding_box(_vertexList, _nVertices);


}


void TriangleMesh::centralize() 
{
    // center object around origin (0,0,0)
    
    // compute mid point of bounding box and translate by its negative vector
    Vector3f mid;
    Vector3f::minus(_boundingBox.getV2(), _boundingBox.getV1(), mid);
    mid.scale(0.5f);
    Vector3f::plus(_boundingBox.getV1(), mid, mid);
    mid.scale(-1.0f);

    this->translate(mid);

    // compute new bounding box
    _boundingBox = compute_bounding_box(_vertexList, _nVertices);
}



Rectangle3f compute_bounding_box(float* vertexList, unsigned int nVertices)
{
    //#undef max 

    float max = 1000.0f;
    
    float mina[3] = {max, max, max};
    float maxa[3] = {-max, -max, -max};


    for (unsigned int i = 0; i < nVertices; i++) {
        if (mina[0] > vertexList[3*i + 0]) {
            mina[0] = vertexList[3*i + 0];
        }

        if (mina[1] > vertexList[3*i + 1]) {
            mina[1] = vertexList[3*i + 1];
        }

        if (mina[2] > vertexList[3*i + 2]) {
            mina[2] = vertexList[3*i + 2];
        }

        if (maxa[0] < vertexList[3*i + 0]) {
            maxa[0] = vertexList[3*i + 0];
        }

        if (maxa[1] < vertexList[3*i + 1]) {
            maxa[1] = vertexList[3*i + 1];
        }

        if (maxa[2] < vertexList[3*i + 2]) {
            maxa[2] = vertexList[3*i + 2];
        }
    }

    Vector3f v1(mina[0], mina[1], mina[2]);
    Vector3f v2(maxa[0], maxa[1], maxa[2]);

    Rectangle3f rect(v1, v2);

    //rect.dump();

    return rect;
}

unsigned int TriangleMesh::getNumFaces() const 
{
    return _nFaces;
}

unsigned int TriangleMesh::getNumVertices() const
{
    return _nVertices;
}

const unsigned int* TriangleMesh::getFaceList() const
{
    return _faceList;
}

const float* TriangleMesh::getVertexList() const
{
    return _vertexList;
}