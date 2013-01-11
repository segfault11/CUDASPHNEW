#include <thrust\sort.h>
#include <thrust\device_ptr.h>
#include <thrust\for_each.h>
#include <thrust\iterator\zip_iterator.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "particle_simulation.h"
#include "util.h"
#include "cgtk\include\clock.h"
#include "boundary_map.h"
#include <thrust\scan.h>

using namespace std;

//-----------------------------------------------------------------------------
//  DEVICE CODE 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  global device variables 
//-----------------------------------------------------------------------------
__constant__ SimulationParameters gSimParamsDev;

texture<float, cudaTextureType1D, cudaReadModeElementType> 
    gParticleVertexData;
texture<float, cudaTextureType1D, cudaReadModeElementType> 
    gParticleSimulationData;
texture<int, cudaTextureType1D, cudaReadModeElementType> 
    gCellStartList;
texture<int, cudaTextureType1D, cudaReadModeElementType> 
    gCellEndList;
texture<int, cudaTextureType1D, cudaReadModeElementType> 
    gSortedParticleIdList;
texture<int, cudaTextureType1D, cudaReadModeElementType> 
    gParticleHashList;

//  information about boundary handling
__constant__ float gBoundaryOrigin[3];
__constant__ float gDx;
__constant__ unsigned int gnBoundarySamples[3];
__constant__ float gRestDistance;

texture<float, cudaTextureType1D, cudaReadModeElementType> 
    gNodeTable;
texture<unsigned int, cudaTextureType1D, cudaReadModeElementType> 
    gIndexMap;

//-----------------------------------------------------------------------------
//  declaration of aux. functions (device) 
//-----------------------------------------------------------------------------
__device__ inline int3 compute_grid_coordinate (float3 pos, float d);
__device__ inline int compute_hash_from_grid_coordinate (int i, int j, int k);
__device__ inline float compute_distance (float3 a, float3 b);
__device__ inline float norm (const float3& a);
__device__ inline void normalize (float3& a);
__device__ inline float dot_product (const float3& a, const float3& b);
__device__ float compute_particle_density_cell (const float3 &pos, 
	float* pParticleList, int* pParticleIdList, int start, int end);
__device__ inline void compute_viscosity_pressure_forces_and_ifsurf_cell
    (const float3& xi, float rhoi, float pi, const float3& vi,
    float* particleVertexData, float* particleSimulationData, 
    int* particleIdList, int start, int end, 
    float3* force, float3* colGra, float* colLapl,
    float3* sumPosNeighbor, float* nNeighbors);

//-----------------------------------------------------------------------------
// CUDA Kernel definitions 
//-----------------------------------------------------------------------------
__global__ void compute_particle_hash (float* particleVertexData, 
    int* particleIdList, int* particleHashList, unsigned int numParticles) 
{
    int idx = blockIdx.x*blockDim.x+threadIdx.x;

    if (idx >= numParticles)
    {
        return;
    }

    // calculate corresponding gridpoint
    int x = (int)((tex1Dfetch(gParticleVertexData, idx*VD_NUM_ELEMENTS + VD_POS_X) - 
		gSimParamsDev.gridOrigin[0])/gSimParamsDev.gridSpacing);
    int y = (int)((tex1Dfetch(gParticleVertexData, idx*VD_NUM_ELEMENTS + VD_POS_Y) - 
		gSimParamsDev.gridOrigin[1])/gSimParamsDev.gridSpacing);
    int z = (int)((tex1Dfetch(gParticleVertexData, idx*VD_NUM_ELEMENTS + VD_POS_Z) - 
		gSimParamsDev.gridOrigin[2])/gSimParamsDev.gridSpacing);

    // wrap outer particles to grid
    // TODO: modulo operation using "&" is faster, requires grid dims of 
	// power of two
    x = x % gSimParamsDev.gridDim[0];
    y = y % gSimParamsDev.gridDim[1];
    z = z % gSimParamsDev.gridDim[2];
    
    // calculate hash, i.e. grid cell id
    int hash = gSimParamsDev.gridDim[0]*(gSimParamsDev.gridDim[1]*z + y) + x;

    particleIdList[idx] = idx;
    particleHashList[idx] = hash;
}
//-----------------------------------------------------------------------------
__global__ void compute_sub_particle_hash (float* particleVertexData, 
    int* particleIdList, int* particleHashList, unsigned int numParticles) 
{
    int idx = blockIdx.x*blockDim.x+threadIdx.x;

    if (idx >= numParticles)
    {
        return;
    }

    int id = particleIdList[idx];

    // calculate corresponding gridpoint
    int x = (int)((particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_X] - 
        gSimParamsDev.gridOrigin[0])/gSimParamsDev.gridSpacingSubParticles);
    int y = (int)((particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Y] - 
		gSimParamsDev.gridOrigin[1])/gSimParamsDev.gridSpacingSubParticles);
    int z = (int)((particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Z] - 
		gSimParamsDev.gridOrigin[2])/gSimParamsDev.gridSpacingSubParticles);

    // wrap outer particles to grid
    // TODO: modulo operation using "&" is faster, requires grid dims of 
	// power of two
    x = x % gSimParamsDev.gridDimSubParticles[0];
    y = y % gSimParamsDev.gridDimSubParticles[1];
    z = z % gSimParamsDev.gridDimSubParticles[2];
    
    // calculate hash, i.e. grid cell id
    int hash = gSimParamsDev.gridDimSubParticles[0]*
        (gSimParamsDev.gridDimSubParticles[1]*z + y) + x;

    particleHashList[idx] = hash;
}
//-----------------------------------------------------------------------------
__global__ void compute_cell_start_end (int* particleHashList, 
	int* cellStartList,  int* cellEndList, unsigned int numParticles)
{
    extern __shared__ int sharedHash[];
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    int hash;

    if (idx < numParticles) 
    {
        hash = particleHashList[idx];
        sharedHash[threadIdx.x + 1] = hash;
        
        if (idx > 0 && threadIdx.x == 0) 
        {
            sharedHash[0] = particleHashList[idx - 1];
        }
    }

    __syncthreads();

    if (idx < numParticles) 
    {
        if (idx == 0 || hash != sharedHash[threadIdx.x])
        {
            cellStartList[hash] = idx;
        
            if (idx > 0) 
            {
                cellEndList[sharedHash[threadIdx.x]] = idx;
            }
        }

        if (idx == numParticles - 1)
        {
            cellEndList[hash] = idx + 1;
        }
    }
}
//-----------------------------------------------------------------------------
//  Compute density and pressure for each particle 
__global__ void compute_particle_density_pressure (float* particleVertexData, 
	float* particleSimulationData, int* particleIdList, int* cellStartList, 
    int* cellEndList) 
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if(idx >= gSimParamsDev.numParticles) {
        return;
    }

	int id = particleIdList[idx];

	float density = 0.0f;
    float pressure;
    float3 pos;

    // get particles position form vertex data
    pos.x = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_X);
    pos.y = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_Y);
    pos.z = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_Z);

    int3 c0 = compute_grid_coordinate(pos, -gSimParamsDev.compactSupport);
    int3 c1 = compute_grid_coordinate(pos, gSimParamsDev.compactSupport);

    int hash;
    int start;
    int end;

    for(int k = c0.z; k <= c1.z; k++) 
    {
        for(int j = c0.y; j <= c1.y; j++) 
        {
            for(int i = c0.x; i <= c1.x; i++)
            {
                hash = compute_hash_from_grid_coordinate(i, j, k);
                start = tex1Dfetch(gCellStartList, hash);
                end = tex1Dfetch(gCellEndList, hash);
                density += compute_particle_density_cell(pos, 
                    particleVertexData, particleIdList, start, end);
            }
        }
    }
    
    density *= gSimParamsDev.particleMass;
    pressure = gSimParamsDev.gasStiffness*(density - 
        gSimParamsDev.restDensity);

    // set density and pressure
    particleSimulationData[id*SD_NUM_ELEMENTS + SD_DENSITY] = density;
    particleSimulationData[id*SD_NUM_ELEMENTS + SD_PRESSURE] = pressure; 
}
//-----------------------------------------------------------------------------
__global__ void compute_particle_acceleration_ifsurf 
    (float* particleVertexData, float* particleSimulationData, 
    int* particleIdList, int* cellStartList,
    int* cellEndList, int* isSurfaceParticle)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }
    
    int id = tex1Dfetch(gSortedParticleIdList, idx);

    float density  = tex1Dfetch(gParticleSimulationData, id*SD_NUM_ELEMENTS
        + SD_DENSITY);
    float pressure = tex1Dfetch(gParticleSimulationData, id*SD_NUM_ELEMENTS
        + SD_PRESSURE);
    float tenCoeff = gSimParamsDev.tensionCoefficient;

    float3 pos;
    pos.x = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_X);
    pos.y = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_Y);
    pos.z = tex1Dfetch(gParticleVertexData, id*VD_NUM_ELEMENTS + VD_POS_Z);

    float3 vel;
    vel.x = tex1Dfetch(gParticleSimulationData, 
        id*SD_NUM_ELEMENTS + SD_VEL0_X);
    vel.y = tex1Dfetch(gParticleSimulationData, 
        id*SD_NUM_ELEMENTS + SD_VEL0_Y);
    vel.z = tex1Dfetch(gParticleSimulationData, 
        id*SD_NUM_ELEMENTS + SD_VEL0_Z);

    int3 c0 = compute_grid_coordinate(pos, -gSimParamsDev.compactSupport);
    int3 c1 = compute_grid_coordinate(pos, gSimParamsDev.compactSupport);

    float3 force;
    force.x = 0.0f;
    force.y = 0.0f;
    force.z = 0.0f;

    float3 colGra;
    colGra.x = 0.0f;
    colGra.y = 0.0f;
    colGra.z = 0.0f;

    // [sumPosNeigbor] and [nNeigbors] are used to computed the center of mass
    // of the neighborhood of this particle (this also includes the particle
    // itself
    float3 sumPosNeighbor;
    sumPosNeighbor.x = pos.x;
    sumPosNeighbor.x = pos.y;
    sumPosNeighbor.x = pos.z;

    float nNeighbors = 1.0f;

    float colLapl;
    float colGraNorm;
    float grav = gSimParamsDev.gravity;

    int hash;
    int start;
    int end;

    // compute viscosity and pressure forces
    for(int k = c0.z; k <= c1.z; k++)
    {
        for(int j = c0.y; j <= c1.y; j++)
        {
            for(int i = c0.x; i <= c1.x; i++)
            {
                hash  = compute_hash_from_grid_coordinate(i, j, k);
                start = tex1Dfetch(gCellStartList, hash);
                end = tex1Dfetch(gCellEndList, hash);
                compute_viscosity_pressure_forces_and_ifsurf_cell(pos, density, 
                    pressure, vel, particleVertexData, particleSimulationData,
                    particleIdList, start, end, &force, &colGra, &colLapl,
                    &sumPosNeighbor, &nNeighbors);
            }
        }
    }

    // surface tension
    colGraNorm = sqrtf(colGra.x*colGra.x + colGra.y*colGra.y 
        + colGra.z*colGra.z);

    float fCoeff = tenCoeff*colLapl/colGraNorm;

    if(colGraNorm > gSimParamsDev.normThresh) 
    {
        force.x -= fCoeff*colGra.x;
        force.y -= fCoeff*colGra.y;
        force.z -= fCoeff*colGra.z;
    }



    // compute contribution of boundary to the pressure force
    /*unsigned int i, j, k;
    
    i = (unsigned int)((pos.x - gBoundaryOrigin[0])/gDx);
    j = (unsigned int)((pos.y - gBoundaryOrigin[1])/gDx);
    k = (unsigned int)((pos.z - gBoundaryOrigin[2])/gDx);

    unsigned int idx2 = i + gnBoundarySamples[0]*(j + gnBoundarySamples[1]*k);
    unsigned int nodeIdx = tex1Dfetch(gIndexMap, idx2);
    float dist = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_DISTANCE);
    
    float3 bNorm;

    bNorm.x = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_X);
    bNorm.y = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_Y);
    bNorm.z = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_Z);
    
    float3 boundaryForce;
    float bCoeff;

    bCoeff = gSimParamsDev.particleMass*(gRestDistance - dist)/
        (gSimParamsDev.timeStep*gSimParamsDev.timeStep);

    boundaryForce.x = bCoeff*bNorm.x;
    boundaryForce.y = bCoeff*bNorm.y;
    boundaryForce.z = bCoeff*bNorm.z;*/

    
    // store the actual acceleration
    particleSimulationData[id*SD_NUM_ELEMENTS + SD_ACC_X] = force.x/density;  
    particleSimulationData[id*SD_NUM_ELEMENTS + SD_ACC_Y] = force.y/density
        - grav;  
    particleSimulationData[id*SD_NUM_ELEMENTS + SD_ACC_Z] = force.z/density;  

    // find out if particle is a surface particle
    /*int isSurface = 0;

    float3 dCenterMass;
    dCenterMass.x = pos.x - sumPosNeighbor.x/nNeighbors;
    dCenterMass.y = pos.y - sumPosNeighbor.y/nNeighbors;
    dCenterMass.z = pos.z - sumPosNeighbor.z/nNeighbors;

    float dCenterMassNormSq = dCenterMass.x*dCenterMass.x +
        dCenterMass.y*dCenterMass.y + dCenterMass.z*dCenterMass.z;

    if (colGraNorm > gSimParamsDev.normThresh) {
        isSurfaceParticle[id] = 1;
    } else {
        isSurfaceParticle[id] = 0;
    }*/
}
//----------------------------------------------------------------------------
/*__global__ void compute_sub_particle_acceleration
    (float* particleVertexData, float* particleSimulationData, 
    int* particleIdList, int* cellStartList,
    int* cellEndList, int* isSurfaceParticle)
{

}*/
//-----------------------------------------------------------------------------
__global__ void integrate_euler (float* particleVertexData, 
    float* particleSimulationData)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }

    unsigned int idVert = idx*VD_NUM_ELEMENTS;
    unsigned int idSim = idx*SD_NUM_ELEMENTS;
    float dt = gSimParamsDev.timeStep;

    particleSimulationData[idSim + SD_VEL0_X] += 
        dt*particleSimulationData[idSim + SD_ACC_X];
    particleSimulationData[idSim + SD_VEL0_Y] += 
        dt*particleSimulationData[idSim + SD_ACC_Y];
    particleSimulationData[idSim + SD_VEL0_Z] += 
        dt*particleSimulationData[idSim + SD_ACC_Z];

    particleVertexData[idVert + VD_POS_X] += 
        dt*particleSimulationData[idSim + SD_VEL0_X];
    particleVertexData[idVert + VD_POS_Y] += 
        dt*particleSimulationData[idSim + SD_VEL0_Y];
    particleVertexData[idVert + VD_POS_Z] += 
        dt*particleSimulationData[idSim + SD_VEL0_Z];   
}
//-----------------------------------------------------------------------------
__global__ void integrate_sub_particles_euler (float* subParticleVertexData, 
    float* subParticleSimulationData, int* subParticleIds, 
    unsigned int nSubParticles)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= nSubParticles)
    {
        return;
    }

    int id = subParticleIds[idx];
    unsigned int idVert = id*VD_NUM_ELEMENTS;
    unsigned int idSim = id*SD_NUM_ELEMENTS;
    float dt = gSimParamsDev.timeStep;
    
    
    subParticleSimulationData[idSim + SD_VEL0_X] += 
        dt*subParticleSimulationData[idSim + SD_ACC_X];
    subParticleSimulationData[idSim + SD_VEL0_Y] += 
        dt*subParticleSimulationData[idSim + SD_ACC_Y];
    subParticleSimulationData[idSim + SD_VEL0_Z] += 
        dt*subParticleSimulationData[idSim + SD_ACC_Z];
    
    subParticleVertexData[idVert + VD_POS_X] += 
        dt*subParticleSimulationData[idSim + SD_VEL0_X];
    subParticleVertexData[idVert + VD_POS_Y] += 
        dt*subParticleSimulationData[idSim + SD_VEL0_Y];
    subParticleVertexData[idVert + VD_POS_Z] += 
        dt*subParticleSimulationData[idSim + SD_VEL0_Z];
}
//-----------------------------------------------------------------------------
__global__ void shift_state (char* particleState)
{
        unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

        if (idx >= gSimParamsDev.numParticles)
        {
            return;        
        }

        particleState[idx] = (particleState[idx] << 2);
}
//-----------------------------------------------------------------------------
__global__ void collision_handling (float* particleVertexData, 
    float* particleSimulationData)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }

    unsigned int idVert = idx*VD_NUM_ELEMENTS;
    unsigned int idSim = idx*SD_NUM_ELEMENTS;

    float3 pos;
    float3 vel;

    pos.x = tex1Dfetch(gParticleVertexData, idVert + VD_POS_X);
    pos.y = tex1Dfetch(gParticleVertexData, idVert + VD_POS_Y);
    pos.z = tex1Dfetch(gParticleVertexData, idVert + VD_POS_Z);

    vel.x = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_X);
    vel.y = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_Y);
    vel.z = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_Z);

    float3 local;
    float3 diff;
    float3 nrm;

    float dist;
    float depth;

    // compute "distance" to box, if positive the particle
    // is outside the box.

    // compute local position of the particle to the box
    local.x = pos.x - gSimParamsDev.boxCen[0];
    local.y = pos.y - gSimParamsDev.boxCen[1];
    local.z = pos.z - gSimParamsDev.boxCen[2];

    // project local pos to the upper right quadrand and
    // compute difference to the boxDim vec
    diff.x = abs(local.x) - gSimParamsDev.boxDim[0];
    diff.y = abs(local.y) - gSimParamsDev.boxDim[1];
    diff.z = abs(local.z) - gSimParamsDev.boxDim[2];

    dist = max(diff.x, diff.y);
    dist = max(dist, diff.z);
    
    // if the particle lies outside the box, the collision must be handled
    float3 contact;
    
    if (dist > 0.0f)
    {
        // contact point in "box space"
        contact.x = min(gSimParamsDev.boxDim[0], 
            max(-gSimParamsDev.boxDim[0], local.x));
        contact.y = min(gSimParamsDev.boxDim[1],
            max(-gSimParamsDev.boxDim[1], local.y));
        contact.z = min(gSimParamsDev.boxDim[2],
            max(-gSimParamsDev.boxDim[2], local.z));

        // translate to worldspace
        contact.x += gSimParamsDev.boxCen[0];
        contact.y += gSimParamsDev.boxCen[1];
        contact.z += gSimParamsDev.boxCen[2];

        // compute penetration depth
        depth = compute_distance(contact, pos);

        // compute normal
        nrm.x = pos.x - contact.x;
        nrm.y = pos.y - contact.y;
        nrm.z = pos.z - contact.z;
        normalize(nrm);

        float velNorm = norm(vel);
        float dp     = dot_product(nrm, vel);
        float coeff  = (1 + gSimParamsDev.restitution*depth/
            (gSimParamsDev.timeStep*velNorm))*dp;

        vel.x -= coeff*nrm.x;
        vel.y -= coeff*nrm.y;
        vel.z -= coeff*nrm.z;

        particleVertexData[idVert + VD_POS_X] = contact.x;
        particleVertexData[idVert + VD_POS_Y] = contact.y;
        particleVertexData[idVert + VD_POS_Z] = contact.z;

        particleSimulationData[idSim + SD_VEL0_X] = vel.x;
        particleSimulationData[idSim + SD_VEL0_Y] = vel.y;
        particleSimulationData[idSim + SD_VEL0_Z] = vel.z;
    }
}
//-----------------------------------------------------------------------------
__global__ void update_particle_state (float* particleVertexData, 
    char* particleState, int* particleIdList, int* cellStartList, 
    int* cellEndList)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }

    unsigned int id = particleIdList[idx];
    float3 pos;
    float3 xj;
    float3 r;
    float rn;
    pos.x = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_X];

    if (pos.x >= 0.2f && pos.x <= 0.5f)
    {
        particleState[id] |= 1;

        // distribute information to neigbors
        pos.y = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Y];
        pos.z = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Z];

        int3 c0 = compute_grid_coordinate(pos, -gSimParamsDev.compactSupport);
        int3 c1 = compute_grid_coordinate(pos, gSimParamsDev.compactSupport);
        
        int hash;
        int start;
        int end;
        
        for(int k = c0.z; k <= c1.z; k++)
        {
            for(int j = c0.y; j <= c1.y; j++)
            {
                for(int i = c0.x; i <= c1.x; i++)
                {
                    hash  = compute_hash_from_grid_coordinate(i, j, k);
                    start = cellStartList[hash];
                    end = cellEndList[hash];
                    
                    for (int u = start; u < end; u++) 
                    {
                        int v = particleIdList[u]; 
                        xj.x = particleVertexData[v*VD_NUM_ELEMENTS + VD_POS_X];
                        xj.y = particleVertexData[v*VD_NUM_ELEMENTS + VD_POS_Y];                         
                        xj.z = particleVertexData[v*VD_NUM_ELEMENTS + VD_POS_Z];
                        r.x = pos.x - xj.x;
                        r.y = pos.y - xj.y;
                        r.z = pos.z - xj.z;
                        rn = r.x*r.x + r.y*r.y + r.z*r.z;

                        if (rn <= gSimParamsDev.compactSupport*
                            gSimParamsDev.compactSupport)
                        {
                            particleState[v] |= 2;
                        }
                    }
                }
            }
        }
        
        // why is this neccessary?
        particleState[id] |= 1;
    }
}
//-----------------------------------------------------------------------------
// sets the sub particle vertex & simulation data for each particle, that has
// changed its state from "default" to "boundary" or "split".
// This kernel is called for particles that were split and for boundary 
// particles. [numParticles] refers to the total number of those particles.
__global__ void initialize_sub_particles (float* subParticleVertexData, 
    float* subParticleSimulationData, int* particleIds, 
    float* particleVertexData, float* particleSimulationData, 
    char* particleStates, unsigned int numParticles)
{
#define SQRT3INV 0.577350269 // = 1/sqrt{3}

    // directions to seed the new sub particles
    const float directions[] = {
                                    SQRT3INV,  SQRT3INV,  SQRT3INV, 
                                    SQRT3INV,  SQRT3INV, -SQRT3INV,
                                    SQRT3INV, -SQRT3INV,  SQRT3INV,
                                    SQRT3INV, -SQRT3INV, -SQRT3INV,
                                   -SQRT3INV,  SQRT3INV,  SQRT3INV,
                                   -SQRT3INV,  SQRT3INV, -SQRT3INV,
                                   -SQRT3INV, -SQRT3INV,  SQRT3INV,
                                   -SQRT3INV, -SQRT3INV, -SQRT3INV
                               };

    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= numParticles)
    {
        return;
    }

    unsigned int id = particleIds[idx];
    char state = particleStates[id];

    // if the base particle was "default" previously ...
    if ((state & 12) == 0)
    {
        float density = particleSimulationData[id*SD_NUM_ELEMENTS + SD_DENSITY];
        float radicand = 3.0f*gSimParamsDev.particleMass/(4.0f*M_PI*density);
        float radius = pow(radicand, 1.0f/3.0f);
        float3 pos;
        pos.x = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_X];
        pos.y = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Y];
        pos.z = particleVertexData[id*VD_NUM_ELEMENTS + VD_POS_Z];

        // ... initialize initial position and velocity of the corresponding 
        // sub particles
        for (unsigned int i = 0; i < 8; i++) 
        {
            // update velocity
            int index = (8*id + i)*SD_NUM_ELEMENTS;
            subParticleSimulationData[index + SD_VEL0_X] = 0.0f;
                particleSimulationData[id*SD_NUM_ELEMENTS + SD_VEL0_X];
            subParticleSimulationData[index + SD_VEL0_Y] = 0.1f;
                particleSimulationData[id*SD_NUM_ELEMENTS + SD_VEL0_Y];
            subParticleSimulationData[index + SD_VEL0_Z] = 0.0f;
                particleSimulationData[id*SD_NUM_ELEMENTS + SD_VEL0_Z];
           
            // update position
            index = (8*id + i)*VD_NUM_ELEMENTS;
            subParticleVertexData[index + VD_POS_X] = pos.x + 
                directions[3*i + 0]*radius;
            subParticleVertexData[index + VD_POS_Y] = pos.y + 
                directions[3*i + 1]*radius;
            subParticleVertexData[index + VD_POS_Z] = pos.z + 
                directions[3*i + 2]*radius;
        }
    }

#undef SQRT3INV
}
//-----------------------------------------------------------------------------
__global__ void check_split_boundary_default (char* particleState, 
        int* isSplit, int* isBoundary, int* isDefault)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    char state = particleState[idx] & 3;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }

    if (state == 0)
    {
        isDefault[idx] = 1;
    }
    else if (state == 2)
    {
        isBoundary[idx] = 1;
    }
    else
    {
        isSplit[idx] = 1;
    }
}
//-----------------------------------------------------------------------------
__global__ void collect_ids (int* subParticleIdList, int* splitParticleIdList,
    int* boundaryParticleIdList, int* defaultParticleIdList,
    int* isSplit, int* isBoundary, int* isDefault, int* splitPrefixSum, 
    int* boundaryPrefixSum, int* defaultPrefixSum, 
    unsigned int numParticlesSplit)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

    if (idx >= gSimParamsDev.numParticles)
    {
        return;
    }

    if (isSplit[idx] == 1)
    {
        int splitPreSum = splitPrefixSum[idx];

        for (unsigned int i = 0; i < 8; i++)
        {
            subParticleIdList[8*splitPreSum + i] = 8*idx + i;
        }

        splitParticleIdList[splitPreSum] = idx; 
    }
    else if (isBoundary[idx] == 1)
    {
        int boundaryPreSum = boundaryPrefixSum[idx];

        for (unsigned int i = 0; i < 8; i++)
        {
            subParticleIdList[8*(numParticlesSplit + boundaryPreSum) + i] =  
                8*idx + i;
        }

        boundaryParticleIdList[boundaryPreSum] = idx;
    }
    else
    {
        defaultParticleIdList[defaultPrefixSum[idx]] = idx;
    }
}
//__global__ void collision_handling(float* particleVertexData, 
//    float* particleSimulationData)
//{
//    /*unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
//
//    if (idx >= gSimParamsDev.numParticles)
//    {
//        return;
//    }
//
//    unsigned int idVert = idx*VD_NUM_ELEMENTS;
//    unsigned int idSim = idx*SD_NUM_ELEMENTS;
//
//    float3 pos;
//    float3 vel;
//
//    pos.x = tex1Dfetch(gParticleVertexData, idVert + VD_POS_X);
//    pos.y = tex1Dfetch(gParticleVertexData, idVert + VD_POS_Y);
//    pos.z = tex1Dfetch(gParticleVertexData, idVert + VD_POS_Z);
//
//    vel.x = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_X);
//    vel.y = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_Y);
//    vel.z = tex1Dfetch(gParticleSimulationData, idSim + SD_VEL0_Z);
//
//    //
//    unsigned int i,j,k;
//    i = (unsigned int)((pos.x - gBoundaryOrigin[0])/gDx);
//    j = (unsigned int)((pos.y - gBoundaryOrigin[1])/gDx);
//    k = (unsigned int)((pos.z - gBoundaryOrigin[2])/gDx);
//    unsigned int idx2 = i + gnBoundarySamples[0]*(j + gnBoundarySamples[1]*k);
//    unsigned int nodeIdx = tex1Dfetch(gIndexMap, idx2);
//    float dist = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_DISTANCE);
//    
//    float3 bNorm;
//
//    bNorm.x = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_X);
//    bNorm.y = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_Y);
//    bNorm.z = tex1Dfetch(gNodeTable, NC_NUM_ELEMENTS*nodeIdx + NC_NORMAL_Z);
//
//    if (bNorm.y != 0.0f)
//    {
//        particleVertexData[idVert + VD_POS_X] -= gSimParamsDev.timeStep*vel.x;
//        particleVertexData[idVert + VD_POS_Y] -= gSimParamsDev.timeStep*vel.y;
//        particleVertexData[idVert + VD_POS_Z] -= gSimParamsDev.timeStep*vel.z;
//    }*/
//}

//-----------------------------------------------------------------------------
// definition of aux. functions (device) 
//-----------------------------------------------------------------------------
__device__ inline int3 compute_grid_coordinate(float3 pos, float d)
{
    int3 gridCoord;

    gridCoord.x = (unsigned int)((pos.x + d - gSimParamsDev.gridOrigin[0])/
        gSimParamsDev.gridSpacing);
    gridCoord.y = (unsigned int)((pos.y + d - gSimParamsDev.gridOrigin[1])/
        gSimParamsDev.gridSpacing);
    gridCoord.z = (unsigned int)((pos.z + d - gSimParamsDev.gridOrigin[2])/
        gSimParamsDev.gridSpacing);

    gridCoord.x = gridCoord.x%gSimParamsDev.gridDim[0];
    gridCoord.y = gridCoord.y%gSimParamsDev.gridDim[1];
    gridCoord.z = gridCoord.z%gSimParamsDev.gridDim[2];

    gridCoord.x = min(max(gridCoord.x, 0),gSimParamsDev.gridDim[0] - 1);
    gridCoord.y = min(max(gridCoord.y, 0),gSimParamsDev.gridDim[1] - 1);
    gridCoord.z = min(max(gridCoord.z, 0),gSimParamsDev.gridDim[2] - 1);

    return gridCoord;
}
//-----------------------------------------------------------------------------
__device__ inline int compute_hash_from_grid_coordinate(int i, int j, int k)
{
    return gSimParamsDev.gridDim[0]*(gSimParamsDev.gridDim[1]*k + j) + i;
}
//-----------------------------------------------------------------------------
__device__ inline float norm(const float3& a)
{
    return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
//-----------------------------------------------------------------------------
__device__ inline void normalize (float3& a)
{
    float norm = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
    a.x /= norm;
    a.y /= norm;
    a.z /= norm;
}
//-----------------------------------------------------------------------------
//  Computes the Euclidean distance between two points.
__device__ inline float compute_distance (float3 a, float3 b)
{
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) 
        + (a.z-b.z)*(a.z-b.z));
}
//-----------------------------------------------------------------------------
__device__ inline float dot_product (const float3& a, const float3& b)  
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
//-----------------------------------------------------------------------------
//  Computes the contribution of neighborparticles of one particular grid cell
//  to the density of the particle at position [pos].
__device__ float compute_particle_density_cell (const float3 &pos, 
	float* particleVertexData, int* particleIdList, int start, int end)
{
    int particleIndex; // index of the neighbor of the particle
    float density = 0.0f;
    float3 p; // neighbor particle's position
    float h = gSimParamsDev.compactSupport;
    float r;
    float d;

    for (int i = start; i < end; i++) 
    {
        particleIndex = tex1Dfetch(gSortedParticleIdList, i);

        // compute position of the neighbor
        p.x = tex1Dfetch(gParticleVertexData, particleIndex*VD_NUM_ELEMENTS
            + VD_POS_X);
        p.y = tex1Dfetch(gParticleVertexData, particleIndex*VD_NUM_ELEMENTS
            + VD_POS_Y);
        p.z = tex1Dfetch(gParticleVertexData, particleIndex*VD_NUM_ELEMENTS
            + VD_POS_Z);

        r = compute_distance(p, pos);
        
        // TODO: evaluating r*r <= h*h might save taking the sqrt in 
        // compute_distance proc. 
        if (r <= h) 
        {
            d = h*h - r*r;
            density += gSimParamsDev.poly6*d*d*d;
        }
    }

    return density;
}
//-----------------------------------------------------------------------------
__device__ inline void compute_viscosity_pressure_forces_and_ifsurf_cell
    (const float3& xi, float rhoi, float pi, const float3& vi,
    float* particleVertexData, float* particleSimulationData, 
    int* particleIdList, int start, int end, 
    float3* force, float3* colGra, float* colLapl,
    float3* sumPosNeighbor, float* nNeighbors)
{
    int j;      // neighbor index in particle list
    float3 xj;  // neighbor particle's position
    float3 vj;  // neighbor particle's velocity
    float rhoj; // neighbor density
    float pj;   // neighbor pressure
    float3 r;   // xi - xj
    float rn;   // ||xi - xj||
    float h = gSimParamsDev.compactSupport; // effective radius
    float grad  = gSimParamsDev.gradSpiky;
    float lapl  = gSimParamsDev.laplVisc;
    float grad2 = gSimParamsDev.gradPoly6;
    float lapl2 = gSimParamsDev.laplPoly6;

    float pressure; // pressure term in the kernel approx
    float rhoi2 = rhoi*rhoi;                    
    float m = gSimParamsDev.particleMass;
    float mu = gSimParamsDev.dynamicViscosity;

    float d; // helper value to avoid arithmetic operations

    for (int i = start; i < end; i++) 
    {
        // get neighbor index from particle list
        j = tex1Dfetch(gSortedParticleIdList, i); 

        // get neighbor particle information
        xj.x = tex1Dfetch(gParticleVertexData, j*VD_NUM_ELEMENTS + VD_POS_X);
        xj.y = tex1Dfetch(gParticleVertexData, j*VD_NUM_ELEMENTS + VD_POS_Y);
        xj.z = tex1Dfetch(gParticleVertexData, j*VD_NUM_ELEMENTS + VD_POS_Z);
        vj.x = tex1Dfetch(gParticleSimulationData, j*SD_NUM_ELEMENTS
            + SD_VEL0_X);
        vj.y = tex1Dfetch(gParticleSimulationData, j*SD_NUM_ELEMENTS
            + SD_VEL0_Y);
        vj.z = tex1Dfetch(gParticleSimulationData, j*SD_NUM_ELEMENTS
            + SD_VEL0_Z);
        rhoj = tex1Dfetch(gParticleSimulationData, j*SD_NUM_ELEMENTS
            + SD_DENSITY);
        pj   = tex1Dfetch(gParticleSimulationData, j*SD_NUM_ELEMENTS
            + SD_PRESSURE);

        r.x = xi.x - xj.x;
        r.y = xi.y - xj.y;
        r.z = xi.z - xj.z;

        rn = norm(r);
        
        // TODO: * masse koennte ausgeklammert werden um multiplikationen
        //         zu sparen.
        //       * generell kann der pressure term in hinblick auf rhoi und
        //         pi vereinfacht werden.
        //       * visc force: mu koennte ausgeklammert werden etc.
        //       * zwei float3's fuer beide kraefte koennten genutzt werden
        //         um die terme zu vereinfachen.
        pressure = rhoi*(pi/rhoi2 + pj/(rhoj*rhoj))*m;

        if (rn <= h && rn > 0.0f)
        {
            // compute pressure force
            d = (h-rn)*(h-rn);

            force->x -= pressure*grad*d/rn*r.x;
            force->y -= pressure*grad*d/rn*r.y;
            force->z -= pressure*grad*d/rn*r.z;
        
            // compute viscosity force
            d = (h-rn);

            force->x += mu*(vj.x-vi.x)*m/rhoj*lapl*d;
            force->y += mu*(vj.y-vi.y)*m/rhoj*lapl*d;
            force->z += mu*(vj.z-vi.z)*m/rhoj*lapl*d;

            // compute color gradient
            d = (h*h-rn*rn)*(h*h-rn*rn);

            colGra->x += m/rhoj*grad2*d*r.x;
            colGra->y += m/rhoj*grad2*d*r.y;
            colGra->z += m/rhoj*grad2*d*r.z;

            // compute color laplacian
            d = (h*h - rn*rn)*(3.0f*h*h - 7.0f*rn*rn);

            *colLapl += m/rhoj*lapl2*d;

            //
            sumPosNeighbor->x += xj.x;
            sumPosNeighbor->y += xj.y;
            sumPosNeighbor->z += xj.z;
            *nNeighbors += 1.0f;
        }
    }
}

//-----------------------------------------------------------------------------
//  HOST CODE
//-----------------------------------------------------------------------------

#define EMPTY_CELL 0xFFFFFFFF

//-----------------------------------------------------------------------------
//  forward declaration of aux. functions
//-----------------------------------------------------------------------------
void create_particle_box (float sx, float sy, float sz, float d, 
    unsigned int numParticles, float** particleVD, float** particleSD,
    unsigned int* numParticlesCreated);
void set_simulation_domain (float xs, float ys, float zs, float xe,
    float ye, float ze, float gridSpacing, float gridSpacingSubParticles,
    SimulationParameters* parameters);
void compute_particle_kernel_invocation_information 
    (unsigned int& nThreadsBlock, unsigned int& nBlocks, 
    unsigned int numParticles);


//-----------------------------------------------------------------------------
//  Definition of ParticleSimulation class 
//-----------------------------------------------------------------------------
/* Set everything to NULL/0
*/
ParticleSimulation::ParticleSimulation (): mParticleVertexData(NULL), 
    mParticleSimulationData(NULL), mParticleVertexDataDevPtr(NULL),
    mParticleSimulationDataDevPtr(NULL), mParticleIdsDevPtr(NULL),
    mParticleHashListDevPtr(NULL), mCellStartListDevPtr(NULL), 
    mCellEndListDevPtr(NULL), mIsSurfaceParticleDevPtr(NULL), mParticleVertexDataVbo(0),
    mNumBlocks(0), mThreadsPerBlock(0), mNumSubParticles(0)
{
    memset(&mParameters, 0, sizeof(SimulationParameters));
}
//-----------------------------------------------------------------------------
ParticleSimulation::~ParticleSimulation() 
{
    this->freeAll();
}
//-----------------------------------------------------------------------------
ParticleSimulation* ParticleSimulation::Example01 () 
{
    // create a particle simulation 
    ParticleSimulation* sim = new ParticleSimulation();

    // create box (cube) of particles
    create_particle_box(-0.65f, -0.45f, -0.25f, 0.5f, 40000, 
        &sim->mParticleVertexData, &sim->mParticleSimulationData,
        &sim->mParameters.numParticles);

    if (sim->mParticleVertexData == NULL || 
        sim->mParticleSimulationData == NULL) 
    {
        THROW_EXCEPTION("Could not allocate memory for particles (Host).");
    }

    // set sph simulation related parameters
    sim->mParameters.kernelParticles = 20;
    sim->mParameters.restDensity = 998.648f;
    sim->mParameters.particleMass = sim->mParameters.restDensity*0.5f*0.5f*0.5f/
        static_cast<float>(sim->mParameters.numParticles);
    sim->mParameters.gasStiffness = 3.0f;
    sim->mParameters.dynamicViscosity = 3.0f;
    sim->mParameters.gravity = 9.81f;
    sim->mParameters.tensionCoefficient = 0.0728f;
    sim->mParameters.normThresh = 15.065f;

    // compute the kernel radius
    float h = powf((3.0f*0.5f*0.5f*0.5f*sim->mParameters.kernelParticles)/
        (4.0f*M_PI*sim->mParameters.numParticles), 1.0f/3.0f);

    sim->mParameters.compactSupport =  h;
    sim->mParameters.poly6 =  315.0f/(64.0f*M_PI*h*h*h*h*h*h*h*h*h);
    sim->mParameters.gradPoly6 = -945.0f/(32.0f*M_PI*h*h*h*h*h*h*h*h*h);
    sim->mParameters.laplPoly6 = -945.0f/(32.0f*M_PI*h*h*h*h*h*h*h*h*h);
    sim->mParameters.gradSpiky = -45.0f/(M_PI*h*h*h*h*h*h);
    sim->mParameters.laplVisc  =  45.0f/(M_PI*h*h*h*h*h*h);
    sim->mParameters.timeStep  = 0.003;
    
    set_simulation_domain(-2.5, -2.5, -2.5, 2.5, 2.5, 2.5, h, h/2.0f,
        &sim->mParameters);

    // set fluid volume
    sim->mParameters.fluidVolume = 0.5f*0.5f*0.5f; 

    // set parameters for boundary handling
    sim->mParameters.restitution = 0.0f;
    sim->mParameters.boxCen[0] = 0.0f;
    sim->mParameters.boxCen[1] = 0.0f;
    sim->mParameters.boxCen[2] = 0.0f;
    sim->mParameters.boxDim[0] = 0.7f;    
    sim->mParameters.boxDim[1] = 0.5f;    
    sim->mParameters.boxDim[2] = 0.3f;    

    // set parameters for new boundary handling
    sim->_boundaryMapFileName = std::string("icosphere.txt");

    // set parameters for surface extraction
    sim->mParameters.cmDistanceThresh = 0.5f;
    sim->mParameters.nPartTresh = 20.0f;
    sim->_leftI = 0.0f;
    sim->_rightI = 1.0f;
    //printf("h %")
    return sim;
}
//-----------------------------------------------------------------------------
int* ParticleSimulation::CreateIsParticleSurfaceList
    (const ParticleSimulation* sim)
{
    int* isSurfaceParticleList = new int[sim->mParameters.numParticles];
    
    CUDA_SAFE_CALL( cudaMemcpy(isSurfaceParticleList, 
        sim->mIsSurfaceParticleDevPtr,
        sizeof(int)*sim->mParameters.numParticles, 
        cudaMemcpyDeviceToHost) );

    int extr = 0;
    for (unsigned int i = 0; i < sim->mParameters.numParticles; i++) {
       extr += isSurfaceParticleList[i];
    }

    printf("%d of %d extracted\n", extr, sim->mParameters.numParticles);
    

    return isSurfaceParticleList;
}
//-----------------------------------------------------------------------------
void ParticleSimulation::FreeIsParticleSurfaceList 
    (int** isSurfaceParticleList)
{
    if (*isSurfaceParticleList == NULL) 
    {
        return;
    }

    delete[] *isSurfaceParticleList;
    *isSurfaceParticleList = NULL;
}
//-----------------------------------------------------------------------------
void ParticleSimulation::freeAll () 
{
    // free host memory
    saveDeleteArray<float>(&mParticleVertexData);
    saveDeleteArray<float>(&mParticleSimulationData);
    
    // free device memory

    // free cuda memory
    cudaSafeFree<float>(&mParticleSimulationDataDevPtr);
    cudaSafeFree<int>(&mParticleIdsDevPtr);
    cudaSafeFree<int>(&mParticleHashListDevPtr);
    cudaSafeFree<int>(&mCellStartListDevPtr);
    cudaSafeFree<int>(&mCellEndListDevPtr);
    cudaSafeFree<int>(&mIsSurfaceParticleDevPtr);
    
    // free OpenGL vertex buffer object
    if (mParticleVertexDataVbo != 0) 
    {
        CUDA_SAFE_CALL( cudaGraphicsUnregisterResource(mGraphicsResources[0]) );
        CUDA_SAFE_CALL( cudaGraphicsUnregisterResource(mGraphicsResources[1]) );
        //cudaGLUnregisterBufferObject(mParticleVertexDataVbo); // <- deprecated
        glDeleteBuffers(1, &mParticleVertexDataVbo);
        glDeleteBuffers(1, &mSubParticleVertexDataVbo);
        mParticleVertexDataVbo = 0;
        mSubParticleVertexDataVbo = 0;
    }
}
//-----------------------------------------------------------------------------
void ParticleSimulation::Init () 
{
    //
    // free device memory, if previously allocated 
    //

    // free cuda memory
    cudaSafeFree<float>(&mParticleSimulationDataDevPtr);
    cudaSafeFree<int>(&mParticleIdsDevPtr);
    cudaSafeFree<int>(&mParticleHashListDevPtr);
    cudaSafeFree<int>(&mCellStartListDevPtr);
    cudaSafeFree<int>(&mCellEndListDevPtr);
    cudaSafeFree<int>(&mIsSurfaceParticleDevPtr);
    
    // free OpenGL vertex buffer object
    if (mParticleVertexDataVbo != 0) 
    {
        CUDA_SAFE_CALL( cudaGraphicsUnregisterResource(mGraphicsResources[0]) );
        CUDA_SAFE_CALL( cudaGraphicsUnregisterResource(mGraphicsResources[1]) );
        glDeleteBuffers(1, &mParticleVertexDataVbo);
        glDeleteBuffers(1, &mSubParticleVertexDataVbo);
        mParticleVertexDataVbo = 0;
        mSubParticleVertexDataVbo = 0;
    }

    //
    // allocate cuda device memory for storing the particles' vertex and
    // simulation data.
    // Vertex data is allocated on device using OpenGL, as it is stored
    // in an vertex buffer object, which is used for rendering later.
    //

    // Simulation data is allocated through cuda.
    CUDA_SAFE_CALL( cudaMalloc(&mParticleSimulationDataDevPtr, 
        mParameters.numParticles*sizeof(float)*SD_NUM_ELEMENTS) );

    // copy initial host data to device
    CUDA_SAFE_CALL( cudaMemcpy(mParticleSimulationDataDevPtr, 
        mParticleSimulationData, 
        mParameters.numParticles*sizeof(float)*SD_NUM_ELEMENTS,
        cudaMemcpyHostToDevice) );
    
    // Vertex data is allocated through a vertex buffer object
    // the vbo is then registered to be used with CUDA
    glGenBuffers(1, &mParticleVertexDataVbo);
    glBindBuffer(GL_ARRAY_BUFFER, mParticleVertexDataVbo);
    glBufferData(GL_ARRAY_BUFFER, 
        mParameters.numParticles*VD_NUM_ELEMENTS*sizeof(float),
          mParticleVertexData, GL_DYNAMIC_COPY);
    CUDA_SAFE_CALL( cudaGraphicsGLRegisterBuffer(&mGraphicsResources[0], 
        mParticleVertexDataVbo, cudaGraphicsMapFlagsNone) );
    //cudaGLRegisterBufferObject(mParticleVertexDataVbo); // <- is deprecated
    
    //
    // alloc & Init additional aux. arrays for nearest neighbor search
    //
    const int* dim = mParameters.gridDim; 
    unsigned int size = dim[0]*dim[1]*dim[2]*sizeof(int);

    CUDA_SAFE_CALL( cudaMalloc(&mCellStartListDevPtr, size) );
    CUDA_SAFE_CALL( cudaMalloc(&mCellEndListDevPtr, size) );

    // set each cell to be empty
    CUDA_SAFE_CALL( cudaMemset(mCellStartListDevPtr, EMPTY_CELL, size) );
    CUDA_SAFE_CALL( cudaMemset(mCellEndListDevPtr, EMPTY_CELL, size) );
     
    CUDA_SAFE_CALL( cudaMalloc(&mParticleIdsDevPtr, 
        mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&mParticleHashListDevPtr, 
        mParameters.numParticles*sizeof(int)) );

    // alloc dev memory for surface particle extraction
    CUDA_SAFE_CALL( cudaMalloc(&mIsSurfaceParticleDevPtr, 
        mParameters.numParticles*sizeof(int)) );

    this->allocateMemoryTwoScale();

    // set up textures, for faster memory look-ups through caching
    // NOTE: VertexData needs to be mapped to get a valid device pointer, 
    //       as it is initial not allocated through CUDA's malloc
    cudaChannelFormatDesc descf = cudaCreateChannelDesc(32, 0, 0, 0,
		cudaChannelFormatKindFloat);
    cudaChannelFormatDesc desci = cudaCreateChannelDesc(32, 0, 0, 0,
		cudaChannelFormatKindSigned);
    cudaChannelFormatDesc descu = cudaCreateChannelDesc(32, 0, 0, 0,
		cudaChannelFormatKindUnsigned);

    CUDA_SAFE_CALL ( cudaBindTexture(0, gParticleSimulationData, 
        mParticleSimulationDataDevPtr, descf, 
        sizeof(float)*SD_NUM_ELEMENTS*mParameters.numParticles) );
    this->map();
    CUDA_SAFE_CALL ( cudaBindTexture(0, gParticleVertexData, 
        mParticleVertexDataDevPtr, descf, 
        sizeof(float)*VD_NUM_ELEMENTS*mParameters.numParticles) );
    this->unmap();
    CUDA_SAFE_CALL ( cudaBindTexture(0, gCellStartList, mCellStartListDevPtr, 
        desci, size) );
    CUDA_SAFE_CALL ( cudaBindTexture(0, gCellEndList, mCellEndListDevPtr, 
        desci, size) );
    CUDA_SAFE_CALL ( cudaBindTexture(0, gSortedParticleIdList, mParticleIdsDevPtr, 
        desci, mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL ( cudaBindTexture(0, gParticleHashList, mParticleHashListDevPtr, 
        desci, mParameters.numParticles*sizeof(int)) );

    // set number of CUDA blocks and threads per blocks for each kernel 
    // invocation
    // NOTE:  - chose different values than 256 to try to get more performance
    //        - make threadsPerBlock and blocks function parameters
    compute_particle_kernel_invocation_information(mThreadsPerBlock, mNumBlocks, 
        mParameters.numParticles);
    //mThreadsPerBlock = mParameters.numParticles < 256 ? 
    //    mParameters.numParticles : 256;
    //mNumBlocks = mParameters.numParticles % mThreadsPerBlock == 0 ?
    //    mParameters.numParticles/mThreadsPerBlock : 
    //    mParameters.numParticles/mThreadsPerBlock + 1; 

    
    //
    // Init boundary handling
    //
    /*std::cout << "loading boundary information ... " << std::endl;
    BoundaryMap bmap("icosphere.txt");
    std::cout << "finished loading" << std::endl;
   
    unsigned int nCoords = bmap.GetNumCoordinates();
    unsigned int totalSamples = bmap.GetNumTotalSamples();

    CUDA_SAFE_CALL( cudaMalloc(&_boundaryMapIndexMapDevPtr, 
        totalSamples*sizeof(unsigned int)) );
    
    CUDA_SAFE_CALL( cudaMalloc(&_boundaryMapNodeTableDevPtr, 
        NC_NUM_ELEMENTS*nCoords*sizeof(float)) );


    CUDA_SAFE_CALL( cudaMemcpy(_boundaryMapIndexMapDevPtr, bmap.GetIndexMap(),
        sizeof(unsigned int)*totalSamples, cudaMemcpyHostToDevice) );

    CUDA_SAFE_CALL( cudaMemcpy(_boundaryMapNodeTableDevPtr, bmap.GetNodeTable(),
        NC_NUM_ELEMENTS*nCoords*sizeof(float), cudaMemcpyHostToDevice) );

    CUDA_SAFE_CALL ( cudaBindTexture(0, gIndexMap, _boundaryMapIndexMapDevPtr, 
        descu, totalSamples*sizeof(unsigned int)) );

    CUDA_SAFE_CALL ( cudaBindTexture(0, gNodeTable, _boundaryMapNodeTableDevPtr, 
        descf, NC_NUM_ELEMENTS*nCoords*sizeof(float)) );

    unsigned int nSamples[3];
    nSamples[0] = bmap.GetIMax();
    nSamples[1] = bmap.GetJMax();
    nSamples[2] = bmap.GetKMax();

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gnBoundarySamples, nSamples, 
		3*sizeof(unsigned int), 0, cudaMemcpyHostToDevice) );

    float origin[3];
    origin[0] = bmap.GetDomain().getV1().getX();
    origin[1] = bmap.GetDomain().getV1().getY();
    origin[2] = bmap.GetDomain().getV1().getZ();

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gBoundaryOrigin, origin, 
		3*sizeof(float), 0, cudaMemcpyHostToDevice) );

    float dx = bmap.GetDx();

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gDx, &dx, 
		sizeof(float), 0, cudaMemcpyHostToDevice) );

    float restDist = bmap.GetRestDistance();

    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gRestDistance, &restDist, 
		sizeof(float), 0, cudaMemcpyHostToDevice) );*/

}
//-----------------------------------------------------------------------------
// allocates and initializes memory needed for the two scale particle 
// simulation
void ParticleSimulation::allocateMemoryTwoScale ()
{
    // create opengl vbo for storing the vertex information of the 
    // sub particles
    glGenBuffers(1, &mSubParticleVertexDataVbo);
    glBindBuffer(GL_ARRAY_BUFFER, mSubParticleVertexDataVbo);
    glBufferData(GL_ARRAY_BUFFER, 8*mParameters.numParticles*VD_NUM_ELEMENTS*
        sizeof(float), NULL, GL_DYNAMIC_COPY);
    CUDA_SAFE_CALL( cudaGraphicsGLRegisterBuffer(&mGraphicsResources[1], 
        mSubParticleVertexDataVbo, cudaGraphicsMapFlagsNone) );

    // create opengl vbo for storing the ids of the particles in 
    // default state
    glGenBuffers(1, &mParticleIdsDefaultVbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mParticleIdsDefaultVbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mParameters.numParticles*sizeof(int),
        NULL, GL_DYNAMIC_COPY);
    CUDA_SAFE_CALL( cudaGraphicsGLRegisterBuffer(&mGraphicsResources[2], 
        mParticleIdsDefaultVbo, cudaGraphicsMapFlagsNone) );

    // create opengl vbo for storing the ids of the active sub particles 
    glGenBuffers(1, &mSubParticleIdsVbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mSubParticleIdsVbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mParameters.numParticles*sizeof(int)*8,
        NULL, GL_DYNAMIC_COPY);
    CUDA_SAFE_CALL( cudaGraphicsGLRegisterBuffer(&mGraphicsResources[3], 
        mSubParticleIdsVbo, cudaGraphicsMapFlagsNone) );

    CUDA_SAFE_CALL( cudaMalloc(&mSubParticleSimulationDataDevPtr, 
        8*mParameters.numParticles*sizeof(float)*SD_NUM_ELEMENTS) );
    CUDA_SAFE_CALL( cudaMemset(mSubParticleSimulationDataDevPtr, 0,
        8*mParameters.numParticles*sizeof(float)*SD_NUM_ELEMENTS) );
    CUDA_SAFE_CALL( cudaMalloc(&mSubParticleHashsDevPtr, 
        8*mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&mParticleStatesDevPtr, 
               mParameters.numParticles*sizeof(char)) );
    CUDA_SAFE_CALL( cudaMemset(mParticleStatesDevPtr, 0, 
        mParameters.numParticles*sizeof(char)) );
    mParticleStates = new char[mParameters.numParticles];    
    CUDA_SAFE_CALL( cudaMalloc(&_isSplitDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&_isBoundaryDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&_isDefaultDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&_splitPrefixSumDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );    
    CUDA_SAFE_CALL( cudaMalloc(&_boundaryPrefixSumDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );    
    CUDA_SAFE_CALL( cudaMalloc(&_defaultPrefixSumDevPtr, 
        (mParameters.numParticles + 1)*sizeof(int)) );    
    CUDA_SAFE_CALL( cudaMalloc(&mParticleIdsDefaultDevPtr, 
        mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&mParticleIdsBoundaryDevPtr, 
        mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&mParticleIdsSplitDevPtr, 
        mParameters.numParticles*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMalloc(&mSubParticleCellStartIdsDevPtr, 
        sizeof(int)*mParameters.gridDimSubParticles[0]*
        mParameters.gridDimSubParticles[1]*
        mParameters.gridDimSubParticles[2]) );
    CUDA_SAFE_CALL( cudaMalloc(&mSubParticleCellEndIdsDevPtr, 
        sizeof(int)*mParameters.gridDimSubParticles[0]*
        mParameters.gridDimSubParticles[1]*
        mParameters.gridDimSubParticles[2]) );
}
//-----------------------------------------------------------------------------
void ParticleSimulation::Bind () const 
{
    // copy simulation parameters to constant memory on device.
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(gSimParamsDev, (void*)&mParameters, 
        sizeof(SimulationParameters)) );  
}
//-----------------------------------------------------------------------------

void ParticleSimulation::Advance ()
{
    this->map();
    this->computeParticleHash();
    this->computeSubParticleHash();
    this->sortParticleIdsByHash();
    this->sortSubParticleIdsByHash();
    this->computeCellStartEndList();
    this->computeSubParticleCellStartEndList ();
    this->computeDensityPressure();
    this->computeAcceleration();
    this->integrate();
    this->integrateSubParticles();
    this->handleCollisions();
    this->computeParticleState();
    this->collect();
    this->initializeSubParticles();
    this->unmap();
}
//-----------------------------------------------------------------------------
float ParticleSimulation::GetParticleRadius () const
{
    return powf((3.0*mParameters.fluidVolume)/
        (4.0*M_PI*mParameters.numParticles), 1.0f/3.0f);
}
//-----------------------------------------------------------------------------
float ParticleSimulation::GetSubParticleRadius () const
{
    return 0.5f*this->GetParticleRadius();
}
//-----------------------------------------------------------------------------
const char* ParticleSimulation::GetParticleState () const
{
    CUDA_SAFE_CALL( cudaMemcpy(mParticleStates, mParticleStatesDevPtr, 
        sizeof(char)*mParameters.numParticles, cudaMemcpyDeviceToHost) );
    return mParticleStates;
}
//-----------------------------------------------------------------------------
unsigned int ParticleSimulation::GetNumParticles () const
{
    return mParameters.numParticles;
}
//-----------------------------------------------------------------------------
void ParticleSimulation::SetNPartThresh (float dVal)
{
    mParameters.nPartTresh += dVal;
    printf("# particle thresh %f\n", mParameters.nPartTresh);
    this->Bind();
}
//-----------------------------------------------------------------------------
void ParticleSimulation::DecreaseCmDistanceThresh ()
{
    _rightI = mParameters.cmDistanceThresh;
    mParameters.cmDistanceThresh = 0.5f*(_rightI - _leftI);
    printf("cmDistance = %f\n", mParameters.cmDistanceThresh);
    this->Bind();
}
//-----------------------------------------------------------------------------
void ParticleSimulation::IncreaseCmDistanceThresh ()
{
    _leftI = mParameters.cmDistanceThresh;
    mParameters.cmDistanceThresh = 0.5f*(_rightI - _leftI);
    printf("cmDistance = %f\n", mParameters.cmDistanceThresh);
    this->Bind();
}
//-----------------------------------------------------------------------------
GLuint ParticleSimulation::GetGLParticleVertexBufferObject () const
{
    return mParticleVertexDataVbo;
}
//-----------------------------------------------------------------------------
GLuint ParticleSimulation::GetGLParticleIndexVertexBufferObject () const
{
    return mParticleIdsDefaultVbo;
}
//-----------------------------------------------------------------------------
unsigned int ParticleSimulation::GetNumParticlesDefault () const
{
    return mNumParticlesDefault;
}
//----------------------------------------------------------------------------
GLuint ParticleSimulation::GetGLSubParticleVertexBufferObject () const
{
    return mSubParticleVertexDataVbo;
}
//----------------------------------------------------------------------------
GLuint ParticleSimulation::GetGLSubParticleIndexVertexBufferObject () const
{
    return mSubParticleIdsVbo;
}
//----------------------------------------------------------------------------
unsigned int ParticleSimulation::GetNumSubParticles () const
{
    return mNumSubParticles;
}
//----------------------------------------------------------------------------
unsigned int ParticleSimulation::GetNumSubParticlesRegular () const
{
    return mNumParticlesSplit*8;
}
//----------------------------------------------------------------------------
unsigned int ParticleSimulation::GetNumSubParticlesBoundary () const
{
    return mNumParticlesBoundary*8;
}
//-----------------------------------------------------------------------------
// Definition of private methods
//-----------------------------------------------------------------------------
void ParticleSimulation::computeParticleHash () 
{
    compute_particle_hash <<< mNumBlocks, mThreadsPerBlock >>> 
        (mParticleVertexDataDevPtr, mParticleIdsDevPtr, 
        mParticleHashListDevPtr, mParameters.numParticles);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeSubParticleHash ()
{
    if (mNumSubParticles != 0)
    {
        compute_sub_particle_hash
            <<< mNumBlocksSubParticle, mThreadsPerBlockSubParticle >>> 
            (mSubParticleVertexDataDevPtr, mSubParticleIdsDevPtr, 
            mSubParticleHashsDevPtr, mNumSubParticles);
    }
}
//-----------------------------------------------------------------------------
void ParticleSimulation::sortParticleIdsByHash ()
{
   thrust::sort_by_key(thrust::device_ptr<int>(mParticleHashListDevPtr),
        thrust::device_ptr<int>(mParticleHashListDevPtr + 
        mParameters.numParticles), 
        thrust::device_ptr<int>(mParticleIdsDevPtr));
}
//-----------------------------------------------------------------------------
void ParticleSimulation::sortSubParticleIdsByHash ()
{
    thrust::sort_by_key(thrust::device_ptr<int>(mSubParticleHashsDevPtr),
        thrust::device_ptr<int>(mSubParticleHashsDevPtr + 
        mNumSubParticles),
        thrust::device_ptr<int>(mSubParticleIdsDevPtr));
}

//-----------------------------------------------------------------------------
void ParticleSimulation::computeCellStartEndList () 
{
    int* dim = mParameters.gridDim; 
    unsigned int size = dim[0]*dim[1]*dim[2]*sizeof(int);

    cudaMemset(mCellStartListDevPtr, EMPTY_CELL, size);
    cudaMemset(mCellEndListDevPtr, EMPTY_CELL, size);
    
    int sharedMemSize = sizeof(int)*(mThreadsPerBlock + 1);
    compute_cell_start_end <<< mNumBlocks, mThreadsPerBlock,  sharedMemSize>>>  
        (mParticleHashListDevPtr, mCellStartListDevPtr, 
        mCellEndListDevPtr, mParameters.numParticles);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeSubParticleCellStartEndList () 
{
    int* dim = mParameters.gridDimSubParticles; 
    unsigned int size = dim[0]*dim[1]*dim[2]*sizeof(int);

    if (mNumSubParticles == 0)
    {
        return;
    }

    cudaMemset(mSubParticleCellStartIdsDevPtr, EMPTY_CELL, size);
    cudaMemset(mSubParticleCellEndIdsDevPtr, EMPTY_CELL, size);
    
    int sharedMemSize = sizeof(int)*(mThreadsPerBlockSubParticle + 1);
    compute_cell_start_end 
        <<< mNumBlocksSubParticle, mThreadsPerBlockSubParticle, sharedMemSize>>>  
        (mSubParticleHashsDevPtr, mSubParticleCellStartIdsDevPtr, 
        mSubParticleCellEndIdsDevPtr, mNumSubParticles);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeDensityPressure () 
{
    compute_particle_density_pressure <<< mNumBlocks, mThreadsPerBlock >>> 
        (mParticleVertexDataDevPtr,mParticleSimulationDataDevPtr, 
        mParticleIdsDevPtr, mCellStartListDevPtr, mCellEndListDevPtr);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeAcceleration ()
{
    compute_particle_acceleration_ifsurf <<< mNumBlocks, mThreadsPerBlock >>> 
        (mParticleVertexDataDevPtr, mParticleSimulationDataDevPtr, 
        mParticleIdsDevPtr, mCellStartListDevPtr, mCellEndListDevPtr,
        mIsSurfaceParticleDevPtr);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeAccelerationSubParticles ()
{
}
//-----------------------------------------------------------------------------
void ParticleSimulation::integrate ()
{
    integrate_euler <<< mNumBlocks, mThreadsPerBlock >>>
        (mParticleVertexDataDevPtr, mParticleSimulationDataDevPtr);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::integrateSubParticles ()
{
    if (mNumSubParticles != 0)
    {
        integrate_sub_particles_euler 
            <<< mNumBlocksSubParticle, mThreadsPerBlockSubParticle >>>
            (mSubParticleVertexDataDevPtr, mSubParticleSimulationDataDevPtr,
            mSubParticleIdsDevPtr, mNumSubParticles);
    }
}
//-----------------------------------------------------------------------------
void ParticleSimulation::handleCollisions ()
{
    collision_handling <<< mNumBlocks, mThreadsPerBlock >>>
        (mParticleVertexDataDevPtr, mParticleSimulationDataDevPtr);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::computeParticleState ()
{
    shift_state <<< mNumBlocks, mThreadsPerBlock >>> (mParticleStatesDevPtr);
    update_particle_state <<< mNumBlocks, mThreadsPerBlock >>> 
        (mParticleVertexDataDevPtr, mParticleStatesDevPtr, mParticleIdsDevPtr,
        mCellStartListDevPtr, mCellEndListDevPtr);

}
//-----------------------------------------------------------------------------
// Collects all id's of active sub particles, particles with state "split",
// particles with state "boundary", particles with that "default" in their own
// arrays and computes the total number of each particle.
void ParticleSimulation::collect ()
{    
    CUDA_SAFE_CALL( cudaMemset(_isSplitDevPtr, 0, 
        (mParameters.numParticles + 1)*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemset(_isBoundaryDevPtr, 0, 
        (mParameters.numParticles + 1)*sizeof(int)) );
    CUDA_SAFE_CALL( cudaMemset(_isDefaultDevPtr, 0, 
        (mParameters.numParticles + 1)*sizeof(int)) );

    check_split_boundary_default <<< mNumBlocks, mThreadsPerBlock >>> 
        (mParticleStatesDevPtr, _isSplitDevPtr, _isBoundaryDevPtr, 
        _isDefaultDevPtr);

    thrust::exclusive_scan(thrust::device_ptr<int>(_isSplitDevPtr),
        thrust::device_ptr<int>(_isSplitDevPtr + mParameters.numParticles
        + 1), thrust::device_ptr<int>(_splitPrefixSumDevPtr));
    thrust::exclusive_scan(thrust::device_ptr<int>(_isBoundaryDevPtr),
        thrust::device_ptr<int>(_isBoundaryDevPtr + mParameters.numParticles
        + 1), thrust::device_ptr<int>(_boundaryPrefixSumDevPtr));
    thrust::exclusive_scan(thrust::device_ptr<int>(_isDefaultDevPtr),
        thrust::device_ptr<int>(_isDefaultDevPtr + mParameters.numParticles
        + 1), thrust::device_ptr<int>(_defaultPrefixSumDevPtr));

    CUDA_SAFE_CALL( cudaMemcpy(&mNumParticlesSplit,
        &_splitPrefixSumDevPtr[mParameters.numParticles], sizeof(int), 
        cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(&mNumParticlesBoundary,
        &_boundaryPrefixSumDevPtr[mParameters.numParticles], sizeof(int), 
        cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(&mNumParticlesDefault,
        &_defaultPrefixSumDevPtr[mParameters.numParticles], sizeof(int), 
        cudaMemcpyDeviceToHost) );    

    collect_ids <<< mNumBlocks, mThreadsPerBlock >>>
        (mSubParticleIdsDevPtr, mParticleIdsSplitDevPtr,
        mParticleIdsBoundaryDevPtr, mParticleIdsDefaultDevPtr, 
        _isSplitDevPtr, _isBoundaryDevPtr, _isDefaultDevPtr, 
        _splitPrefixSumDevPtr, _boundaryPrefixSumDevPtr, 
        _defaultPrefixSumDevPtr, mNumParticlesSplit);

    mNumSubParticles = 8*(mNumParticlesSplit + mNumParticlesBoundary);

    // compute how many cuda blocks and how many threads a block are needed for
    // split particles, boundary particles, default particles.
    compute_particle_kernel_invocation_information(mThreadsPerBlockSplit, 
        mNumBlocksSplit, mNumParticlesSplit);
    compute_particle_kernel_invocation_information(mThreadsPerBlockBoundary, 
        mNumBlocksBoundary, mNumParticlesBoundary);
    compute_particle_kernel_invocation_information(mThreadsPerBlockDefault, 
        mNumBlocksDefault, mNumParticlesDefault);
    compute_particle_kernel_invocation_information(mThreadsPerBlockSubParticle,
        mNumBlocksSubParticle, mNumSubParticles);
}
//-----------------------------------------------------------------------------
//  initializes new sub particles if a parent particle has changed its state 
//  from "default" to "boundary" or "split"
void ParticleSimulation::initializeSubParticles () 
{
    if (mNumParticlesSplit > 0)
    {
        initialize_sub_particles <<<mNumBlocksSplit, mThreadsPerBlockSplit>>> 
            (mSubParticleVertexDataDevPtr, mSubParticleSimulationDataDevPtr,
            mParticleIdsSplitDevPtr, mParticleVertexDataDevPtr, 
            mParticleSimulationDataDevPtr, mParticleStatesDevPtr, 
            mNumParticlesSplit);    
    }

    if (mNumParticlesBoundary > 0)
    {
        initialize_sub_particles <<<mNumBlocksBoundary, mThreadsPerBlockBoundary>>> 
            (mSubParticleVertexDataDevPtr, mSubParticleSimulationDataDevPtr,
            mParticleIdsBoundaryDevPtr, mParticleVertexDataDevPtr, 
            mParticleSimulationDataDevPtr, mParticleStatesDevPtr, 
            mNumParticlesBoundary);
    }
}
//-----------------------------------------------------------------------------
void ParticleSimulation::map () 
{
    cudaGraphicsMapResources(4, mGraphicsResources);
    size_t nBytes;
    cudaGraphicsResourceGetMappedPointer
        (reinterpret_cast<void**>(&mParticleVertexDataDevPtr), &nBytes,
        mGraphicsResources[0]);
    cudaGraphicsResourceGetMappedPointer
        (reinterpret_cast<void**>(&mSubParticleVertexDataDevPtr), &nBytes,
        mGraphicsResources[1]);
    cudaGraphicsResourceGetMappedPointer
        (reinterpret_cast<void**>(&mParticleIdsDefaultDevPtr), &nBytes,
        mGraphicsResources[2]);
    cudaGraphicsResourceGetMappedPointer
        (reinterpret_cast<void**>(&mSubParticleIdsDevPtr), &nBytes,
        mGraphicsResources[3]);
}

void ParticleSimulation::unmap () 
{
    cudaGraphicsUnmapResources(4, mGraphicsResources);
    //cudaGLUnmapBufferObject(mParticleVertexDataVbo);
}
//-----------------------------------------------------------------------------
void ParticleSimulation::SaveInfoTable (const std::string& filename) 
{
    using namespace std;

    ofstream file;

    file.open(filename);

    int* pIdList = new int[mParameters.numParticles];
    int* pHashList = new int[mParameters.numParticles];

    int cellListSize = mParameters.gridDim[0]*mParameters.gridDim[1]*
        mParameters.gridDim[2];

    int* pCellStartList = new int[cellListSize];
    int* pCellEndList = new int[cellListSize];

    //this->map();

    cudaMemcpy(pHashList, mParticleHashListDevPtr, 
        mParameters.numParticles*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(pIdList, mParticleIdsDevPtr, 
        mParameters.numParticles*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(pCellStartList, mCellStartListDevPtr, 
        cellListSize*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(pCellEndList, mCellEndListDevPtr, 
        cellListSize*sizeof(int), cudaMemcpyDeviceToHost);
    
    file << "Number of particles " << mParameters.numParticles << endl; 
    file << setw(8) << "index" << setw(12) << " id" << setw(12) << 
        " hash" << setw(12) << " start" << setw(12) << " end" << endl;

    for (unsigned int i = 0; i < cellListSize; i++) 
    {
        file << setw(8) << i;

        if(i < mParameters.numParticles)
        {
            file << setw(12) << pIdList[i];
            file << setw(12) << pHashList[i];
        } 
        else 
        {
            file << setw(12) << "";
            file << setw(12) << "";
        }

        if(pCellStartList[i] == EMPTY_CELL) 
        {
            file << setw(12) << "";
        } 
        else
        {
            file << setw(12) << pCellStartList[i];
        }
        
        if(pCellEndList[i] == EMPTY_CELL)
        {
            file << setw(12) << "" << endl;
        } 
        else
        {
            file << setw(12) << pCellEndList[i] << endl;
        }
    }

    delete[] pIdList;
    delete[] pHashList;
    delete[] pCellStartList;
    delete[] pCellEndList;

    file.close();

    //this->unmap();
}
//-----------------------------------------------------------------------------
void ParticleSimulation::SaveParticleInfo (const std::string& filename)
{
    using namespace std;

    this->map();

    ofstream file;
    
    file.open(filename);

    float* particleVertexData = 
        new float[VD_NUM_ELEMENTS*mParameters.numParticles]; 
    float* particleSimulationData = 
        new float[SD_NUM_ELEMENTS*mParameters.numParticles]; 

    // copy particle information from device to host
    cudaMemcpy(particleVertexData, mParticleVertexDataDevPtr, 
		VD_NUM_ELEMENTS*mParameters.numParticles*sizeof(float), 
		cudaMemcpyDeviceToHost);
    cudaMemcpy(particleSimulationData, mParticleSimulationDataDevPtr, 
		SD_NUM_ELEMENTS*mParameters.numParticles*sizeof(float), 
		cudaMemcpyDeviceToHost);    
    
    // set max. chars for each column of the table
    int columnWidth = 20;
    
    file << setw(columnWidth) << "Index";
    file << setw(columnWidth) << "X";
    file << setw(columnWidth) << "Y";
    file << setw(columnWidth) << "Z";
    file << setw(columnWidth) << "Density";
    file << setw(columnWidth) << "Pressure";
    file << setw(columnWidth) << "Acc X";
    file << setw(columnWidth) << "Acc Y";
    file << setw(columnWidth) << "Acc Z";
    file << endl;

    for (unsigned int i = 0; i < mParameters.numParticles; i++) 
    {
        file << setw(columnWidth) << i;
        file << setw(columnWidth) 
            << particleVertexData[VD_NUM_ELEMENTS*i + VD_POS_X];
        file << setw(columnWidth) 
            << particleVertexData[VD_NUM_ELEMENTS*i + VD_POS_Y];
        file << setw(columnWidth) 
            << particleVertexData[VD_NUM_ELEMENTS*i + VD_POS_Z];
        file << setw(columnWidth) 
            << particleSimulationData[SD_NUM_ELEMENTS*i + SD_DENSITY];
        file << setw(columnWidth) 
            << particleSimulationData[SD_NUM_ELEMENTS*i + SD_PRESSURE];
        file << setw(columnWidth) 
            << particleSimulationData[SD_NUM_ELEMENTS*i + SD_ACC_X];
        file << setw(columnWidth) 
            << particleSimulationData[SD_NUM_ELEMENTS*i + SD_ACC_Y];
        file << setw(columnWidth) 
            << particleSimulationData[SD_NUM_ELEMENTS*i + SD_ACC_Z];
        // TODO: rest of the params.
        file << endl;
    }


    delete[] particleVertexData;
    delete[] particleSimulationData;

    file.close();

    this->unmap();
}
//-----------------------------------------------------------------------------
/*
unsigned int ParticleSimulation::GetSizeMemoryGPU () const
{
    float size;
    
    size += mNumPar


}*/
//-----------------------------------------------------------------------------
//  definition of aux. functions
//-----------------------------------------------------------------------------
// Creates a set of particles, that are aligned in a cube, given the starting
// point of the box [sx, sy, sz] the length of the cube in each direction [d]
// and the approximate amount of total particles [numParticles].
//
// Returns a pointer to the vertex data of the particles in [particleVD] and
// a pointer to the simulation data of the particles in [particleSD] and the
// actual amount of particles created.
void create_particle_box (float sx, float sy, float sz, float d, 
    unsigned int numParticles, float** particleVD, float** particleSD,
    unsigned int* numParticlesCreated)
{
    // computed number of particles in each direction
    unsigned int num = pow(static_cast<double>(numParticles), 1.0/3.0);
    *numParticlesCreated = num*num*num;

    *particleVD = new float[*numParticlesCreated*VD_NUM_ELEMENTS];
    *particleSD = new float[*numParticlesCreated*SD_NUM_ELEMENTS];

    // check if new failed.
    if ((*particleSD) == NULL || (*particleSD) == NULL)
    {
        *numParticlesCreated = 0;
        return;
    }

    // compute spatial increment
    float dx = d/static_cast<float>(num - 1);

    // seed the particles inside the cube
    
    // set the position of each particle
    unsigned int idx;
    
    for (unsigned int k = 0; k < num; k++) 
    {
		for (unsigned int j = 0; j < num; j++) 
        {
			for (unsigned int i = 0; i < num; i++) 
            {
			    idx = VD_NUM_ELEMENTS*(num*(num*k+j)+i);
                (*particleVD)[idx + VD_POS_X] = sx + i*dx;
                (*particleVD)[idx + VD_POS_Y] = sy + j*dx;
                (*particleVD)[idx + VD_POS_Z] = sz + k*dx;
            }
		}
	}
    
    // set other particles attributes to 0.0f
    memset((*particleSD), 0, 
        sizeof(float)*SD_NUM_ELEMENTS*(*numParticlesCreated));
}
//-----------------------------------------------------------------------------
// Sets the simulation domain in the [parameters], based on a starting point
// [xs, ys, zs] an ending point [xe, ye, ze] and the distance between two
// grid points [gridSpacing].
void set_simulation_domain (float xs, float ys, float zs, float xe,
    float ye, float ze, float gridSpacing, float gridSpacingSubParticles,
    SimulationParameters* parameters)
{
    parameters->gridOrigin[0] = xs;
    parameters->gridOrigin[1] = ys;
    parameters->gridOrigin[2] = zs;
    parameters->gridDim[0] = static_cast<int>((xe - xs)/gridSpacing + 0.5);
    parameters->gridDim[1] = static_cast<int>((ye - ys)/gridSpacing + 0.5);
    parameters->gridDim[2] = static_cast<int>((ze - zs)/gridSpacing + 0.5);
    parameters->gridDimSubParticles[0] = 
        static_cast<int>((xe - xs)/gridSpacingSubParticles + 0.5);
    parameters->gridDimSubParticles[1] = 
        static_cast<int>((ye - ys)/gridSpacingSubParticles + 0.5);
    parameters->gridDimSubParticles[2] = 
        static_cast<int>((ze - zs)/gridSpacingSubParticles + 0.5);
    parameters->gridSpacing = gridSpacing;
    parameters->gridSpacingSubParticles = gridSpacingSubParticles;
}
//-----------------------------------------------------------------------------
void compute_particle_kernel_invocation_information
    (unsigned int& nThreadsBlock, unsigned int& nBlocks, 
    unsigned int numParticles)
{
    if (numParticles == 0)
    {
        nThreadsBlock = 0;
        nBlocks = 0;
        return;
    }

    nThreadsBlock = numParticles > 256 ? 256 : numParticles;
    nBlocks = numParticles % nThreadsBlock == 0 ? numParticles/nThreadsBlock : 
        numParticles/nThreadsBlock + 1;
}