//-----------------------------------------------------------------------------
//  SphInComplexShapes.cpp
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------

#include "SphInComplexShapes.h"
#include <Wm5DistPoint3Rectangle3.h>
#include <Wm5DistPoint3Box3.h>
#include "Wm5Matrix3.h"
#include "portable_pixmap.h"
#include <stdexcept>

//-----------------------------------------------------------------------------
// Macros
//-----------------------------------------------------------------------------
#define M_PI       3.14159265358979323846
//-----------------------------------------------------------------------------
// Declaration of file intern aux. functions. 
//-----------------------------------------------------------------------------
inline bool operator< (const Wm5::Vector3f& s, const Wm5::Vector3f& e);
inline float evaluateDensityKernel (float sqDist, float h);
inline float evaluateViscosityKernel (float sqDist, float h);
inline float computeSquaredDistance(const Wm5::Vector3f& a,
    const Wm5::Vector3f& b);
inline bool isInsideBox (const Wm5::Vector3f& p, const Wm5::Box3f& b);
inline float computeSignedDistancePoint3Box3(const Wm5::Vector3f& pos, 
    const Wm5::Box3f& box);
//-----------------------------------------------------------------------------
// Public class method definitions
//-----------------------------------------------------------------------------
SphInComplexShapes::SphInComplexShapes (Wm5::Vector3f s, Wm5::Vector3f e, 
        float gridSpacing, float restDistance, float compactSupport, 
        float particleMass, float restDensity, float viscosity,
        float particleSpacing)
: mGridStart(s), mGridEnd(e), mGridSpacing(gridSpacing), 
    mRestDistance(restDistance), mCompactSupport(compactSupport),
    mParticleMass(particleMass), mRestDensity(restDensity), 
    mViscosity(viscosity), mParticleSpacing(particleSpacing)
{
    if (!(s < e))
    {
        throw std::runtime_error("In SphInComplexShapes::SphInComplexShapes: \
         invalid domain definition");
    }

    mGridDimensions[0] = static_cast<unsigned int>((e[0] - s[0])/
        mGridSpacing) + 1;
    mGridDimensions[1] = static_cast<unsigned int>((e[1] - s[1])/
        mGridSpacing) + 1;
    mGridDimensions[2] = static_cast<unsigned int>((e[2] - s[2])/
        mGridSpacing) + 1;

    mNumGridNodes = mGridDimensions[0]*mGridDimensions[1]*
        mGridDimensions[2];

    mIndexGrid = new unsigned int[mNumGridNodes];

    mNormals.push_back(Wm5::Vector3f(0.0f, 0.0f, 0.0f));
    mDistances.push_back(restDistance);
    mDensities.push_back(0.0f);
}
//-----------------------------------------------------------------------------
SphInComplexShapes::~SphInComplexShapes ()
{
    delete[] mIndexGrid;
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::SetRectangle (Wm5::Rectangle3f r, Wm5::Vector3f n)
{
    reset();

    // compute "bounding box" of the rectangle
    float d = mCompactSupport + mGridSpacing;
    Wm5::Vector3f dV(d, d, d);
    Wm5::Vector3f start(r.GetMMCorner());
    start -= dV;
    Wm5::Vector3f end(r.GetPPCorner());
    end += dV;

    // compute grid coordinates that need to be altered
    Wm5::Tuple<3, int> cs = computeGridCoordinate(start);
    Wm5::Tuple<3, int> ce = computeGridCoordinate(end);

    // compute distance and normal for each coordinate to be altered
    for (unsigned int k = cs[2]; k <= ce[2]; k++)
    {
        for (unsigned int j = cs[1]; j <= ce[1]; j++)
        {
            for (unsigned int i = cs[0]; i <= ce[0]; i++)
            {
                Wm5::Tuple<3, int> currentCoordinate;
                currentCoordinate[0] = i;
                currentCoordinate[1] = j;
                currentCoordinate[2] = k;
                Wm5::Vector3f currentPosition = 
                    computePositionFromGridCoordinate(currentCoordinate);
                Wm5::DistPoint3Rectangle3f dist(currentPosition, r);
                float distance = dist.Get();
                
                if (distance > mRestDistance)
                {
                    distance = mRestDistance;
                }
                
                // compute normal
                float u = dist.GetRectangleCoordinate(0);
                float v = dist.GetRectangleCoordinate(1);
                Wm5::Vector3f closest = r.Center + r.Axis[0]*u + r.Axis[1]*v;
                Wm5::Vector3f normal = currentPosition - closest;
                normal.Normalize(0.0000001);

                if (n.Dot(normal) >= 0)
                {
                    mDistances.push_back(distance);               
                    mNormals.push_back(normal);
                }
                else
                {
                    mDistances.push_back(-distance);
                    mNormals.push_back(normal*-1.0f);                                   
                }

                // update index grid
                unsigned int index = computeIndexFromGridCoordinate
                    (currentCoordinate);
                mIndexGrid[index] = mDistances.size() - 1;
            }
        }
    }

    // initialize densities to 0.0f
    for (unsigned int i = 1; i < mDistances.size(); i++)
    {
        mDensities.push_back(0.0f);
    }

    // seed ghost particles
    std::vector<Wm5::Vector3f> particlePositions;
    Wm5::Vector3f particlePosition;
    start = r.GetMMCorner();
    unsigned int numStepsU = std::ceil(2*r.Extent[0]/mParticleSpacing);
    unsigned int numStepsV = std::ceil(2*r.Extent[1]/mParticleSpacing);
    unsigned int numStepsNormal = std::ceil(mCompactSupport/mParticleSpacing);

    for (unsigned int k = 0; k <= numStepsNormal; k++)
    {
        for (unsigned int j = 0; j <= numStepsV; j++)
        {
            for (unsigned int i = 0; i <= numStepsU; i++)
            {
                particlePosition = start + (-n*k + r.Axis[1]*j + r.Axis[0]*i)*
                    mParticleSpacing;
                particlePositions.push_back(particlePosition);
            }
        }
    }

    // project each density contribution of the ghost particle to the grid 
    // cells
    for (unsigned int l = 0; l < particlePositions.size(); l++)
    {
        Wm5::Vector3f position = particlePositions[l];
        Wm5::Vector3f delta(mCompactSupport, mCompactSupport, mCompactSupport);
        Wm5::Vector3f start = position - delta;
        Wm5::Vector3f end = position + delta;
        Wm5::Tuple<3, int> cs = computeGridCoordinate(start);
        Wm5::Tuple<3, int> ce = computeGridCoordinate(end);

        for (unsigned int k = cs[2]; k <= ce[2]; k++)
        {
            for (unsigned int j = cs[1]; j <= ce[1]; j++)
            {    
                for (unsigned int i = cs[0]; i <= ce[0]; i++)
                {
                   Wm5::Vector3f nodePosition = 
                        computePositionFromGridCoordinate(i, j, k); 
                    unsigned int index = 
                        computeIndexFromGridCoordinate(i, j, k);
                    unsigned int node = mIndexGrid[index];

                    if (node != 0)
                    {
                        float sqDist = computeSquaredDistance(position, 
                            nodePosition);
                        mDensities[node] += evaluateDensityKernel
                            (sqDist, mCompactSupport);
                    }
                }
            }
        }
    }

    // multiply densities by particle mass
    for (unsigned int i = 1; i < mDistances.size(); i++)
    {
        mDensities[i] *= mParticleMass;
    }
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::SetBox (const Wm5::Box3f& b)
{
    this->reset();
    
    // compute signed distances
    for (unsigned int i = 0; i < mGridDimensions[0]; i++)
    {
        for (unsigned int j = 0; j < mGridDimensions[1]; j++)
        {
            for (unsigned int k = 0; k < mGridDimensions[2]; k++)
            {
                Wm5::Vector3f pos = computePositionFromGridCoordinate(i, j, k);
                float dist = computeSignedDistancePoint3Box3(pos, b);

                if (std::abs(dist) < mCompactSupport)
                {
                    if (std::abs(dist) > mRestDistance)
                    {
                        dist = dist/std::abs(dist)*mRestDistance;
                    }

                    mDistances.push_back(dist);
                    unsigned int idx = computeIndexFromGridCoordinate(i, j, k);
                    mIndexGrid[idx] = mDistances.size() - 1; 
                }
            }
        }    
    }
    
    // seed particles
    float d = mCompactSupport + mParticleSpacing;
    Wm5::Box3f outerBox(b);
    outerBox.Extent[0] += d;
    outerBox.Extent[1] += d;
    outerBox.Extent[2] += d;
    Wm5::Vector3f s(outerBox.Center.X() - outerBox.Extent[0],
        outerBox.Center.Y() - outerBox.Extent[1],
        outerBox.Center.Z() - outerBox.Extent[2]);
    std::vector<Wm5::Vector3f> particlePositions;
    unsigned int iMax = 2*outerBox.Extent[0]/mParticleSpacing;
    unsigned int jMax = 2*outerBox.Extent[1]/mParticleSpacing;
    unsigned int kMax = 2*outerBox.Extent[2]/mParticleSpacing;

    for (unsigned int k = 0; k <= kMax; k++)
    {
        for (unsigned int j = 0; j <= jMax; j++)
        {            
            for (unsigned int i = 0; i <= iMax; i++)
            {
                Wm5::Vector3f pos = s + (outerBox.Axis[0]*i + outerBox.Axis[1]*j
                    + outerBox.Axis[2]*k)*mParticleSpacing;
                if (!isInsideBox(pos, b))
                {
                    particlePositions.push_back(pos);
                }
            }
        }
    }

    // initialize densitiy and viscosity contributions to 0.0f
    for (unsigned int i = 1; i < mDistances.size(); i++)
    {
        mDensities.push_back(0.0f);
        mViscosities.push_back(0.0f);
    }

    // project each density contribution of the ghost particle to the grid 
    // cells
    for (unsigned int l = 0; l < particlePositions.size(); l++)
    {
        Wm5::Vector3f position = particlePositions[l];
        Wm5::Vector3f delta(mCompactSupport, mCompactSupport, mCompactSupport);
        Wm5::Vector3f start = position - delta;
        Wm5::Vector3f end = position + delta;
        Wm5::Tuple<3, int> cs = computeGridCoordinate(start);
        Wm5::Tuple<3, int> ce = computeGridCoordinate(end);

        for (unsigned int k = cs[2]; k <= ce[2]; k++)
        {
            for (unsigned int j = cs[1]; j <= ce[1]; j++)
            {    
                for (unsigned int i = cs[0]; i <= ce[0]; i++)
                {
                   Wm5::Vector3f nodePosition = 
                        computePositionFromGridCoordinate(i, j, k); 
                    unsigned int index = 
                        computeIndexFromGridCoordinate(i, j, k);
                    unsigned int node = mIndexGrid[index];

                    if (node != 0)
                    {
                        float sqDist = computeSquaredDistance(position, 
                            nodePosition);
                        mDensities[node] += evaluateDensityKernel
                            (sqDist, mCompactSupport);
                        mViscosities[node] += evaluateViscosityKernel
                            (sqDist, mCompactSupport);      
                    }    
                }
            }
        }
    }

    // multiply densities by particle mass
    for (unsigned int i = 1; i < mDistances.size(); i++)
    {
        mDensities[i] *= mParticleMass;
        mViscosities[i] *= mViscosity*mParticleMass/mRestDensity; 
    }

}
//-----------------------------------------------------------------------------
const unsigned int* SphInComplexShapes::GetIndexGrid () const
{
    return mIndexGrid;
}
//-----------------------------------------------------------------------------
const float* SphInComplexShapes::GetDistances () const
{
    return mDistances.data();
}
//-----------------------------------------------------------------------------
const float* SphInComplexShapes::GetDensities () const
{
    return mDensities.data();
}
//-----------------------------------------------------------------------------
unsigned int SphInComplexShapes::GetGridDimension (unsigned int i) const
{
    if (i > 2)
    {
        throw std::runtime_error("In SphInComplexShapes::GetGridDimension: \
            i should be 0, 1 or 2");
    }

    return mGridDimensions[i];
}
//-----------------------------------------------------------------------------
unsigned int SphInComplexShapes::GetNumGridNodes () const 
{
    return mNumGridNodes;
}
//-----------------------------------------------------------------------------
const Wm5::Vector3f& SphInComplexShapes::GetGridStart () const
{
    return mGridStart;
}
//-----------------------------------------------------------------------------
unsigned int SphInComplexShapes::GetNumGridNodesTaken () const
{
    return mDistances.size();
}
//-----------------------------------------------------------------------------
float SphInComplexShapes::GetGridSpacing () const
{
    return mGridSpacing;
}
//-----------------------------------------------------------------------------
float SphInComplexShapes::GetRestDistance () const
{
    return mRestDistance;
}
//-----------------------------------------------------------------------------
float* SphInComplexShapes::CreateDistanceTextureData 
    (const SphInComplexShapes& s)
{
    float* texData = new float[s.GetNumGridNodes()];
    
    for (unsigned int k = 0; k < s.mGridDimensions[2]; k++)
    {
        for (unsigned int j = 0; j < s.mGridDimensions[1]; j++)
        {
            for (unsigned int i = 0; i < s.mGridDimensions[0]; i++)
            {
                unsigned int index = s.computeIndexFromGridCoordinate(i, j, k);
                unsigned int nodeIdx = s.mIndexGrid[index];
                texData[index] = s.mDistances[nodeIdx];
            }
        }
    }
    
    return texData;
}
//-----------------------------------------------------------------------------
float* SphInComplexShapes::CreateDensityTextureData 
    (const SphInComplexShapes& s)
{
    float* texData = new float[s.GetNumGridNodes()];
    
    for (unsigned int k = 0; k < s.mGridDimensions[2]; k++)
    {
        for (unsigned int j = 0; j < s.mGridDimensions[1]; j++)
        {
            for (unsigned int i = 0; i < s.mGridDimensions[0]; i++)
            {
                unsigned int index = s.computeIndexFromGridCoordinate(i, j, k);
                unsigned int nodeIdx = s.mIndexGrid[index];
                texData[index] = s.mDensities[nodeIdx];
            }
        }
    }
    
    return texData;
}
//-----------------------------------------------------------------------------
float* SphInComplexShapes::CreateViscosityTextureData 
    (const SphInComplexShapes& s)
{
    float* texData = new float[s.GetNumGridNodes()];
    
    for (unsigned int k = 0; k < s.mGridDimensions[2]; k++)
    {
        for (unsigned int j = 0; j < s.mGridDimensions[1]; j++)
        {
            for (unsigned int i = 0; i < s.mGridDimensions[0]; i++)
            {
                unsigned int index = s.computeIndexFromGridCoordinate(i, j, k);
                unsigned int nodeIdx = s.mIndexGrid[index];
                texData[index] = s.mViscosities[nodeIdx];
            }
        }
    }
    
    return texData;
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::SaveSlicedDistanceMapToPpm 
    (const std::string& filename) const
{
    PortablePixmap ppm(mGridDimensions[0], mGridDimensions[1], 255);
    unsigned int k = mGridDimensions[2]/2;

    for (unsigned int i = 0; i < mGridDimensions[0]; i++)
    {
        for (unsigned int j = 0; j < mGridDimensions[1]; j++)
        {
            Wm5::Tuple<3, int> currentCoordinate;
            currentCoordinate[0] = i;
            currentCoordinate[1] = j;
            currentCoordinate[2] = k;
            unsigned int index = computeIndexFromGridCoordinate
                (currentCoordinate);
            
            float distance = mDistances[mIndexGrid[index]];
            ppm.setJET(i, j, distance/(mRestDistance));
        
        }
    }
    
    ppm.save(filename);
}
//-----------------------------------------------------------------------------
float SphInComplexShapes::ComputeMaxDensity () const
{
    float maxDensity = 0.0f;

    for (unsigned int i = 0; i < mDensities.size(); i++)
    {
        float density = mDensities[i];

        if (density > maxDensity)
        {
            maxDensity = density;
        } 
    }

    return maxDensity;
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::SaveSlicedDensityMapToPpm 
    (const std::string& filename) const
{
    PortablePixmap ppm(mGridDimensions[0], mGridDimensions[1], 255);
    unsigned int k = mGridDimensions[2]/2;
    float maxDensity = 0.0f;

    for (unsigned int i = 0; i < mDensities.size(); i++)
    {
        float density = mDensities[i];

        if (density > maxDensity)
        {
            maxDensity = density;
        } 
    }

    for (unsigned int i = 0; i < mGridDimensions[0]; i++)
    {
        for (unsigned int j = 0; j < mGridDimensions[1]; j++)
        {
            Wm5::Tuple<3, int> currentCoordinate;
            currentCoordinate[0] = i;
            currentCoordinate[1] = j;
            currentCoordinate[2] = k;
            unsigned int index = computeIndexFromGridCoordinate
                (currentCoordinate);

            if (mIndexGrid[index] == 0)
            {
                ppm.setJET(i, j, 0.0f);
            }
            else
            {
                float density = mDensities[mIndexGrid[index]];
                ppm.setJET(i, j, density/maxDensity);
            }
        
        }
    }
    
    ppm.save(filename);
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::SaveSlicedViscosityMapToPpm
    (const std::string& filename) const
{
    PortablePixmap ppm(mGridDimensions[0], mGridDimensions[1], 255);
    unsigned int k = mGridDimensions[2]/2;
    float maxViscosity = 0.0f;

    for (unsigned int i = 0; i < mViscosities.size(); i++)
    {
        float viscosity = mViscosities[i];

        if (viscosity > maxViscosity)
        {
            maxViscosity = viscosity;
        } 
    }

    for (unsigned int i = 0; i < mGridDimensions[0]; i++)
    {
        for (unsigned int j = 0; j < mGridDimensions[1]; j++)
        {
            Wm5::Tuple<3, int> currentCoordinate;
            currentCoordinate[0] = i;
            currentCoordinate[1] = j;
            currentCoordinate[2] = k;
            unsigned int index = computeIndexFromGridCoordinate
                (currentCoordinate);

            if (mIndexGrid[index] == 0)
            {
                ppm.setJET(i, j, 0.0f);
            }
            else
            {
                float viscosity = mViscosities[mIndexGrid[index]];
                ppm.setJET(i, j, viscosity/maxViscosity);
            }
        
        }
    }
    
    ppm.save(filename);
}
//-----------------------------------------------------------------------------
// Private class method definitions
//-----------------------------------------------------------------------------
Wm5::Tuple<3, int> SphInComplexShapes::computeGridCoordinate 
    (const Wm5::Vector3f& p) const
{
    Wm5::Tuple<3, int> coordinate;
    
    coordinate[0] = static_cast<int>((p.X() - mGridStart.X())/mGridSpacing);    
    coordinate[1] = static_cast<int>((p.Y() - mGridStart.Y())/mGridSpacing);    
    coordinate[2] = static_cast<int>((p.Z() - mGridStart.Z())/mGridSpacing);    

    // clamp coordinates if neccessary
    if (coordinate[0] < 0)
    {
        coordinate[0] = 0;
    }

    if (coordinate[0] >= mGridDimensions[0])
    {
        coordinate[0] = mGridDimensions[0] - 1;
    }

    if (coordinate[1] < 0)
    {
        coordinate[1] = 0;
    }

    if (coordinate[1] >= mGridDimensions[1])
    {
        coordinate[1] = mGridDimensions[1] - 1;
    }

    if (coordinate[2] < 0)
    {
        coordinate[2] = 0;
    }

    if (coordinate[2] >= mGridDimensions[2])
    {
        coordinate[2] = mGridDimensions[2] - 1;
    }

    return coordinate;
}
//-----------------------------------------------------------------------------
Wm5::Vector3f SphInComplexShapes::computePositionFromGridCoordinate 
    (const Wm5::Tuple<3, int>& c) const
{
    float pos[3];

    pos[0] = mGridStart.X() + mGridSpacing*c[0];
    pos[1] = mGridStart.Y() + mGridSpacing*c[1];
    pos[2] = mGridStart.Z() + mGridSpacing*c[2];

    return Wm5::Vector3f(pos[0], pos[1], pos[2]);
}
//-----------------------------------------------------------------------------
Wm5::Vector3f SphInComplexShapes::computePositionFromGridCoordinate 
    (int i, int j, int k) const
{
    float pos[3];

    pos[0] = mGridStart.X() + mGridSpacing*i;
    pos[1] = mGridStart.Y() + mGridSpacing*j;
    pos[2] = mGridStart.Z() + mGridSpacing*k;

    return Wm5::Vector3f(pos[0], pos[1], pos[2]);
}
//-----------------------------------------------------------------------------
unsigned int SphInComplexShapes::computeIndexFromGridCoordinate
    (const Wm5::Tuple<3, int>& c) const
{
    return c[0] + mGridDimensions[0]*(c[1] + mGridDimensions[1]*c[2]);
}
//-----------------------------------------------------------------------------
unsigned int SphInComplexShapes::computeIndexFromGridCoordinate
    (int i, int j, int k) const
{
    return i + mGridDimensions[0]*(j + mGridDimensions[1]*k);
}
//-----------------------------------------------------------------------------
void SphInComplexShapes::reset ()
{
    memset(mIndexGrid, 0, sizeof(unsigned int)*mNumGridNodes);
    mNormals.clear();
    mDistances.clear();
    mDensities.clear();
    mViscosities.clear();
    mNormals.push_back(Wm5::Vector3f(0.0f, 0.0f, 0.0f));
    mDistances.push_back(-mRestDistance);
    mDensities.push_back(0.0f);
    mViscosities.push_back(0.0f);
}
//-----------------------------------------------------------------------------
// Definition of file intern aux. functions. 
//-----------------------------------------------------------------------------
bool operator< (const Wm5::Vector3f& s, const Wm5::Vector3f& e)
{
    if (e.X() < s.X() || e.Y() < s.Y() || e.Z() < s.Z())
    {
        return false;
    }

    return true;
}
//-----------------------------------------------------------------------------
float evaluateDensityKernel (float sqDist, float h)
{
    float coeff = 315.0f/(64.0f*M_PI*h*h*h*h*h*h*h*h*h);
    float sqh = h*h;
    float diff = sqh - sqDist;
    
    if (diff < 0.0f)
    {
        return 0.0f;
    }

    return coeff*diff*diff*diff;
}
//-----------------------------------------------------------------------------
float evaluateViscosityKernel (float sqDist, float h)
{
    float dist = std::sqrtf(sqDist);
    
    if (dist > h)
    {
        return 0.0f;
    }

    return 45.0f/(M_PI*h*h*h*h*h*h)*(h - dist);
}
//-----------------------------------------------------------------------------
float computeSquaredDistance (const Wm5::Vector3f& a, const Wm5::Vector3f& b)
{
    Wm5::Vector3f d = a - b;

    return d.SquaredLength();
}
//-----------------------------------------------------------------------------
bool isInsideBox (const Wm5::Vector3f& p, const Wm5::Box3f& b)
{
    Wm5::Matrix3f m(b.Axis[0], b.Axis[1], b.Axis[2], false);
    Wm5::Vector3f t = p - b.Center;
    t = m*t;

    if (std::abs(t.X()) >= b.Extent[0] || std::abs(t.Y()) >= b.Extent[1] ||
        std::abs(t.Z()) >= b.Extent[2])
    {
        return false;
    }

    return true;
 }
//-----------------------------------------------------------------------------
float computeSignedDistancePoint3Box3(const Wm5::Vector3f& pos,
    const Wm5::Box3f& box)
{
    Wm5::Matrix3f m(box.Axis[0], box.Axis[1], box.Axis[2], false);
    Wm5::Vector3f t = pos - box.Center;
    t = m*t;

    float a = std::abs(t.X()) - box.Extent[0];
    float b = std::abs(t.Y()) - box.Extent[1];
    float c = std::abs(t.Z()) - box.Extent[2];

    float distance;

    if (std::abs(a) < std::abs(b))
    {
        distance = a;
    }
    else
    {
        distance = b;
    }

    if (std::abs(distance) > std::abs(c))
    {
        distance = c;
    }

    if (a > 0.0f || b > 0.0f || c > 0.0f)
    {
        distance = Wm5::DistPoint3Box3f(pos, box).Get();
    }

    return distance;
}
//-----------------------------------------------------------------------------
