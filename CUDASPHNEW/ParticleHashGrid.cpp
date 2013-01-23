//-----------------------------------------------------------------------------
//  ParticleHashGrid.cpp
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------

#include "ParticleHashGrid.h"


//-----------------------------------------------------------------------------
//  public class member definitions
//-----------------------------------------------------------------------------
ParticleHashGrid::ParticleHashGrid (const Wm5::Vector3f& start, 
    const Wm5::Vector3f& end, float gridSpacing)
: mStart(start), mEnd(end), mGridSpacing(gridSpacing)
{
    mGridDimensions[0] = static_cast<unsigned int>(end.X() - start.X()/
        mGridSpacing); 
    mGridDimensions[1] = static_cast<unsigned int>(end.Y() - start.Y()/
        mGridSpacing); 
    mGridDimensions[2] = static_cast<unsigned int>(end.Z() - start.Z()/
        mGridSpacing); 

    mNumGridCells = mGridDimensions[0]*mGridDimensions[1]*
        mGridDimensions[2];

    mHashGrid = new std::set<Wm5::Vector3f>[mNumGridCells];
}
//-----------------------------------------------------------------------------
ParticleHashGrid::~ParticleHashGrid ()
{
    for (unsigned int i = 0; i < mNumGridCells; i++)
    {
        mHashGrid[i].clear();
    }

    delete[] mHashGrid;
}
//-----------------------------------------------------------------------------
void ParticleHashGrid::AddParticlePosition (const Wm5::Vector3f& position) 
{
    Wm5::Tuple<3, int> coordinate = computeCoordinate(position);
    unsigned int index = computeIndex(coordinate);
    mHashGrid[index].insert(position);
}
//-----------------------------------------------------------------------------
void ParticleHashGrid::RangeQuery (const std::set<Wm5::Vector3f> result, 
        const Wm5::Vector3f& query, float range) const
{
    Wm5::Vector3f start;
    start.X() = query.X() - range; 
    start.Y() = query.Y() - range; 
    start.Z() = query.Z() - range; 
    
    Wm5::Vector3f end;
    end.X() = query.X() + range; 
    end.Y() = query.Y() + range; 
    end.Z() = query.Z() + range; 

    Wm5::Tuple<3, int> sc =  computeCoordinate(start);
    Wm5::Tuple<3, int> se =  computeCoordinate(end);

    for (int k = sc[2]; k < ec[2]; k++)
    {
        for (int j = sc[1]; j < ec[1]; j++)
        {
            for (int i = sc[0]; i < ec[0]; i++)
            {
            
            }
        }
    }


}
//-----------------------------------------------------------------------------
//  private (inline) class member definitions
//-----------------------------------------------------------------------------
unsigned int ParticleHashGrid::computeIndex 
    (const Wm5::Tuple<3, int>&  coordinates) const
{
    Wm5::Tuple<3, int> c(coordinates);
    
    // clamp as boundary handling
    if (c[0] < 0)
    {
        c[0] = 0;
    }

    if (c[0] >= mGridDimensions[0])
    {
        c[0] = mGridDimensions[0] - 1;
    }

    if (c[1] < 0)
    {
        c[1] = 0;
    }

    if (c[1] >= mGridDimensions[1])
    {
        c[1] = mGridDimensions[1] - 1;
    }

    if (c[2] < 0)
    {
        c[2] = 0;
    }

    if (c[2] >= mGridDimensions[2])
    {
        c[2] = mGridDimensions[2] - 1;
    }

    return coordinates[0] + mGridDimensions[0]*(coordinates[1] + 
        mGridDimensions[1]*coordinates[2]);
}
//-----------------------------------------------------------------------------
Wm5::Tuple<3, int> ParticleHashGrid::computeCoordinate
    (const Wm5::Vector3f& position) const
{
    Wm5::Tuple<3, int> coordinate;

    coordinate[0] = (position.X() - mStart.X())/mGridSpacing;
    coordinate[1] = (position.Y() - mStart.Y())/mGridSpacing;
    coordinate[2] = (position.Z() - mStart.Z())/mGridSpacing;

    return coordinate;
}
//-----------------------------------------------------------------------------
