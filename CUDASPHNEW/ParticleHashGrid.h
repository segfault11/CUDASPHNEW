//-----------------------------------------------------------------------------
//  ParticleHashGrid.h
//  CUDASPHNEW
//
//  Created by Arno in Wolde Lübke on 22.01.13.
//  Copyright (c) 2013. All rights reserved.
//-----------------------------------------------------------------------------
#ifndef _PARTICLEHASHGRID_H
#define _PARTICLEHASHGRID_H

#include <set>
#include <Wm5Tuple.h>
#include <Wm5Vector3.h>

class ParticleHashGrid
{
public:
    ParticleHashGrid (const Wm5::Vector3f& start, const Wm5::Vector3f& end, 
        float gridSpacing);
    ~ParticleHashGrid ();

    void AddParticlePosition (const Wm5::Vector3f& position);
    void RangeQuery (const std::set<Wm5::Vector3f> result, 
        const Wm5::Vector3f& query, float range) const ;

private:
    ParticleHashGrid ();
    ParticleHashGrid (const ParticleHashGrid& orig);
    ParticleHashGrid& operator= (const ParticleHashGrid& orig);

    inline unsigned int computeIndex 
        (const Wm5::Tuple<3, int>& coordinates) const;
    inline Wm5::Tuple<3, int> computeCoordinate
        (const Wm5::Vector3f& position) const;

    Wm5::Vector3f mStart;
    Wm5::Vector3f mEnd;
    float mGridSpacing;

    std::set<Wm5::Vector3f>* mHashGrid;
    unsigned int mGridDimensions[3];
    unsigned int mNumGridCells;
};

#endif // ParticleHashGrid.h