#ifndef _PARTICLE_SIMULATION_H
#define _PARTICLE_SIMULATION_H

#include <Windows.h>
#include <gl\glew.h>
#include <string>
#include <cuda_gl_interop.h>
#include "util.h"
#include "SphInComplexShapes.h"
#include <iostream>
using namespace std;

/* Enumeration of the vertex data of each particle (i.e. the information which 
** used for rendering the particle and is stored in the OpenGL vertex buffer
** object).
*/
enum ParticleVertexDataIdx
{
    VD_POS_X, 
    VD_POS_Y, 
    VD_POS_Z,
    VD_NUM_ELEMENTS
};

/* Enumeration of additional information about the particles which is needed 
** for the physics simulation.
*/
enum ParticleSimulationDataIdx
{
    SD_MASS,
    SD_DENSITY,
    SD_PRESSURE,
    SD_ACC_X,
    SD_ACC_Y,
    SD_ACC_Z,
    SD_VEL0_X,
    SD_VEL0_Y,
    SD_VEL0_Z,
	SD_ENERGY,
    SD_NUM_ELEMENTS
};


enum BoundaryParticleData
{
    BPD_POS_X,
    BPD_POS_Y,
    BPD_POS_Z,
    BPD_DENSITY,
    BPD_PRESSURE
};


struct SimulationParameters
{
    // Definition of the computational domain ( = uniform cartesian grid) 
    // for simulating particles and sub particles.
    float gridOrigin[3];
    float gridSpacing;
    int gridDim[3]; // # of grid cell in each direction
    float gridSpacingSubParticles;
    int gridDimSubParticles[3];

    /* Particle System info */
    unsigned int numParticles;    /* total amount of particles */

    /* SPH info */
    float compactSupport;
    float compactSupportSub;
    int kernelParticles;        /* Average number of particles within kernel 
                                ** radius */

    /* Integration info */
    float timeStep;
    float timeStepSubParticles;

    /* Simulation info */
    float particleMass;         /* Mass of each particle */
    float subParticleMass;
    float gasStiffness;
    float restDensity;
    float dynamicViscosity;
    float gravity;              /* Acceleration due to gravity */
    float tensionCoefficient;
    float normThresh;

    float fluidVolume;

    // kernel coefficients for particle and sub-particle sim
    float poly6;
    float gradPoly6;
    float laplPoly6;
    float gradSpiky;
    float laplVisc;
    float poly6Sub;
    float gradPoly6Sub;
    float laplPoly6Sub;
    float gradSpikySub;
    float laplViscSub;
    

    /* surface extraction threshold parameters */
    float cmDistanceThresh;
    float nPartTresh;

    /* provisional collision handling */
    float boxDim[3];            /* Vector pointing from the box center to the 
                                ** upper right point of the box */
    float boxCen[3];            /* Center of the box in worldspace. */
    float restitution;
};

class ParticleSimulation
{
public:
    ~ParticleSimulation ();

    // Allocates memory to store particle information on the device and copies
    // the initial particle information to it.
    void Init ();

    // Copies the simulation parameters of the particle simulation to the 
    // device.
    // Should be used before if another particle simulation was active before.
    void Bind () const;

    // Executes one physics step/advances the particle system one step in time
    void Advance ();

    // Executes one physics step only on the (regular) sub particles (debug method)
    void AdvanceSubParticles ();

    void AdvanceTwoScale ();


    // Debug Method to check if the 3d textures are working
    void Check3DTextures () const;


    // Access to the opengl vertex buffer object that stores the vertex data
    // (as described in the [ParticleVertexDataIdx] enum) of the particles
    GLuint GetGLParticleVertexBufferObject () const;

    // Access to the opengl vertex buffer that stores the indices of the 
    // currently not split base particles
    GLuint GetGLParticleIndexVertexBufferObject () const;

    // Access to the number of particles that are currently not split
    unsigned int GetNumParticlesDefault () const;
    unsigned int GetNumSubParticles () const;
    unsigned int GetNumSubParticlesRegular () const;
    unsigned int GetNumSubParticlesBoundary () const;
    GLuint GetGLSubParticleVertexBufferObject () const;
    GLuint GetGLSubParticleIndexVertexBufferObject () const;

    float GetParticleRadius () const;
    float GetSubParticleRadius () const;
    const char* GetParticleState () const;
    unsigned int GetNumParticles () const;
    unsigned int GetNumTimesSteps () const;
    void SetNPartThresh (float val);
    void IncreaseCmDistanceThresh ();
    void DecreaseCmDistanceThresh ();
    void SaveInfoTable (const std::string& filename);
    void SaveParticleInfo (const std::string& filename);

    unsigned int GetSizeMemoryGPU () const;

    /** @brief Creates an examples particle simulation.
    ***
    *** In this examples particles for a cube of fluid that is being dropped in 
    *** a rectangular container.
    **/
    static ParticleSimulation* Example01 ();
    
    /** @brief Creates a list, that indicates whether a particle belongs to \
    ***        the surface layer or not.
    ***
    *** Creates a list, that indicates whether a particle of simulation belongs 
    *** to the surface layer or not. The list resides in host memory.
    **/
    static int* CreateIsParticleSurfaceList (const ParticleSimulation* sim);

    /** @brief Frees an isSurfaceParticleList
    **/
    static void FreeIsParticleSurfaceList (int** isSurfaceParticleList);

private:

    //
    // convenience locks
    //
    
    ParticleSimulation ();
    ParticleSimulation (const ParticleSimulation& orig);
    ParticleSimulation& operator= (const ParticleSimulation& orig);
    
    //
    // Private methods
    //

    // Frees all memory allocated by the object.
    void freeAll ();

    // Map and unmap vertex buffer object to CUDA memory space.
    inline void map ();
    inline void unmap ();
    inline void allocateMemoryTwoScale ();
    inline void computeParticleHash ();
    inline void computeSubParticleHash ();
    inline void sortParticleIdsByHash ();
    inline void sortSubParticleIdsByHash ();
    inline void computeCellStartEndList ();
    inline void computeSubParticleCellStartEndList ();
    inline void computeDensityPressure ();
    inline void computeSubParticleDensityPressure(); 
    inline void computeAcceleration ();
    inline void computeSubParticleAcceleration (); 
    inline void projectQuantities ();
    inline void integrate ();
    inline void integrateSubParticles ();
    inline void handleCollisions ();
    inline void handleSubParticleCollisions ();
    inline void extractSurfaceParticles ();
    inline void computeParticleState ();
    inline void collect ();
    inline void initializeSubParticles ();
    inline void setUpSphInComplexShapes ();

    // Member declarations

    // Host information
    float* mParticleVertexData;             
    float* mParticleSimulationData;         
    char* mParticleStates;

    // CUDA/OpenGL interoperation information 
    GLuint mParticleVertexDataVbo;
    GLuint mParticleIdsDefaultVbo;
    GLuint mSubParticleVertexDataVbo;
    GLuint mSubParticleIdsVbo;
    GLuint mSurfaceParticlesVbo;
    cudaGraphicsResource_t mGraphicsResources[4];
    unsigned int mNumSurfaceParticles;

    // linear arrays stored on the device
    float* mParticleVertexDataDevPtr;         
    float* mParticleSimulationDataDevPtr;
    char* mParticleStatesDevPtr;
    int* mParticleHashListDevPtr;       
    int* mCellStartListDevPtr;          
    int* mCellEndListDevPtr;         
    int* mIsSurfaceParticleDevPtr;
    int* mParticleIdsDevPtr;
    int* _isSplitDevPtr;
    int* _isBoundaryDevPtr;
    int* _isDefaultDevPtr;
    int* _splitPrefixSumDevPtr;
    int* _boundaryPrefixSumDevPtr;
    int* _defaultPrefixSumDevPtr;
    int* mSubParticleIdsDevPtr;
    int* mSubParticleSortedIdsDevPtr;
    int* mParticleIdsDefaultDevPtr;
    int* mParticleIdsBoundaryDevPtr;
    int* mParticleIdsSplitDevPtr;
    float* mSubParticleVertexDataDevPtr; 
    float* mSubParticleSimulationDataDevPtr;
    int* mSubParticleHashsDevPtr;
    int* mSubParticleCellStartIdsDevPtr;          
    int* mSubParticleCellEndIdsDevPtr;     

    // kernel invocation information
    unsigned int mThreadsPerBlock;          
    unsigned int mThreadsPerBlockSplit;
    unsigned int mThreadsPerBlockBoundary;
    unsigned int mThreadsPerBlockDefault;
    unsigned int mThreadsPerBlockSubParticle;
    unsigned int mThreadsPerBlockSubParticleBoundary;
    unsigned int mThreadsPerBlockSubParticleRegular;
    unsigned int mNumBlocks;                   
    unsigned int mNumBlocksSplit;
    unsigned int mNumBlocksBoundary;
    unsigned int mNumBlocksDefault;
    unsigned int mNumBlocksSubParticle;
    unsigned int mNumBlocksSubParticleBoundary;
    unsigned int mNumBlocksSubParticleRegular;

    // numbers of different particle types
    int mNumParticlesSplit;
    int mNumParticlesBoundary;
    int mNumParticlesDefault;
    int mNumSubParticles;

    /* Host and device information */
    SimulationParameters mParameters;       /* simulation parameters */

    /* boundary handling information */
    std::string _boundaryMapFileName;
    float* _boundaryMapNodeTableDevPtr;
    unsigned int* _boundaryMapIndexMapDevPtr;

    /* */
    float _leftI;
    float _rightI;



    unsigned int mNumTimeSteps;
    CudaTimer mTimer;


    float* mDeviceMemory;

    // boundary handling
    float* mBoundaryParticleDataDevPtr;
    unsigned int mBoundaryIndicesDevPtr;
    int* mBoundaryParticleCellStartListDevPtr;          
    int* mBoundaryParticleCellEndListDevPtr;

    // boundary handling / sph in complex shapes
    SphInComplexShapes* mBoundaryHandling;
    cudaArray* mBoundaryDistances;
    cudaArray* mBoundaryDensities;
    cudaArray* mBoundaryViscosities;
};

#endif /*include guard of: particle_simulation.h */