//-----------------------------------------------------------------------------
//  SphInComplexShapes.h
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------
#ifndef _SPHINCOMPLEXSHAPES_H
#define _SPHINCOMPLEXSHAPES_H

#include <vector>
#include <string>
#include <Wm5Vector3.h>
#include <Wm5Rectangle3.h>
#include <Wm5Tuple.h>
#include <Wm5Box3.h>

class SphInComplexShapes
{
public:
    SphInComplexShapes (Wm5::Vector3f s, Wm5::Vector3f e, float gridSpacing, 
        float restDistance, float compactSupport, float particleMass, 
        float restDensity, float viscosity, float particleSpacing);
    ~SphInComplexShapes ();

    //! Sets a plane as an obstacle to the boundary constraints.
    //! \param p Plane.
    //! \param n Normal of the plane. Pointing in negative direction.
    void SetRectangle (Wm5::Rectangle3f r, Wm5::Vector3f n);
    void SetRectangle (Wm5::Vector3f s, Wm5::Vector3f e, Wm5::Vector3f n);
    void SetBox (const Wm5::Box3f& b);

    // class access
    const unsigned int* GetIndexGrid () const;
    const float* GetDistances () const;
    const float* GetDensities () const;
    unsigned int GetGridDimension (unsigned int i) const;
    unsigned int GetNumGridNodes () const;
    const Wm5::Vector3f& GetGridStart () const;
    unsigned int GetNumGridNodesTaken () const;
    float GetGridSpacing () const;
    float GetRestDistance () const;

    float ComputeMaxDensity () const;

    static float* CreateDistanceTextureData (const SphInComplexShapes& s);
    static float* CreateDensityTextureData (const SphInComplexShapes& s);
    static float* CreateViscosityTextureData (const SphInComplexShapes& s);

    void SaveSlicedDistanceMapToPpm (const std::string& filename) const;
    void SaveSlicedDensityMapToPpm (const std::string& filename) const;
    void SaveSlicedViscosityMapToPpm (const std::string& filename) const;

private:
    SphInComplexShapes ();
    SphInComplexShapes (const SphInComplexShapes& orig);
    SphInComplexShapes& operator= (const SphInComplexShapes& orig);

    // Grid related auxiliary functions
    inline Wm5::Tuple<3, int> computeGridCoordinate 
        (const Wm5::Vector3f& p) const;
    inline Wm5::Vector3f computePositionFromGridCoordinate 
        (const Wm5::Tuple<3, int>& c) const;
    inline Wm5::Vector3f computePositionFromGridCoordinate 
        (int i, int j, int k) const;
    inline unsigned int computeIndexFromGridCoordinate 
        (const Wm5::Tuple<3, int>& c) const;
    inline unsigned int computeIndexFromGridCoordinate
        (int i, int j, int k) const;

    // reset the boundary index grid and the node data lists
    inline void reset ();

    float mParticleMass;
    float mParticleSpacing;
    float mGridSpacing;
    float mRestDistance;
    float mCompactSupport;
    float mRestDensity;
    float mViscosity;

    Wm5::Vector3f mGridStart;
    Wm5::Vector3f mGridEnd;
    unsigned int mNumGridNodes;
    unsigned int mGridDimensions[3];
    unsigned int* mIndexGrid;
    
    // node data lists
    std::vector<Wm5::Vector3f> mNormals;
    std::vector<float> mDistances;    
    std::vector<float> mDensities;
    std::vector<float> mViscosities;
};

#endif // SPHBoundaryConstraints.h