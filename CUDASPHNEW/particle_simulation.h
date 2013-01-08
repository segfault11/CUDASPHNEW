#ifndef _PARTICLE_SIMULATION_H
#define _PARTICLE_SIMULATION_H

#include <Windows.h>
#include <gl\glew.h>
#include <string>
#include <cuda_gl_interop.h>
#include "util.h"

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

struct SimulationParameters
{
    /* Grid info */
    float gridOrigin[3];
    float gridSpacing;
    int gridDim[3];             /* Amount of _GRIDCELLS_(!) in each direction */

    /* Particle System info */
    unsigned int nParticles;    /* total amount of particles */

    /* SPH info */
    float compactSupport;
    int kernelParticles;        /* Average number of particles within kernel 
                                ** radius */

    /* Integration info */
    float timeStep;

    /* Simulation info */
    float particleMass;         /* Mass of each particle */
    float gasStiffness;
    float restDensity;
    float dynamicViscosity;
    float gravity;              /* Acceleration due to gravity */
    float tensionCoefficient;
    float normThresh;

    float fluidVolume;

    /* Kernel coefficients */
    float poly6;
    float gradPoly6;
    float laplPoly6;
    float gradSpiky;
    float laplVisc;

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

    void Advance ();
    GLuint GetGLVertexBufferObject () const;
    float GetParticleRadius () const;
    float GetSubParticleRadius () const;
    const char* GetParticleState () const;
    unsigned int GetNumParticles () const;
    void SetNPartThresh (float val);
    void IncreaseCmDistanceThresh ();
    void DecreaseCmDistanceThresh ();
    void SaveInfoTable (const std::string& filename);
    void SaveParticleInfo (const std::string& filename);


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
    ParticleSimulation& operator = (const ParticleSimulation& orig);
    
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
    inline void sortParticleIdsByHash ();
    inline void computeCellStartEndList ();
    inline void computeDensityPressure ();
    inline void computeAcceleration ();
    inline void integrate ();
    inline void handleCollisions ();
    inline void extractSurfaceParticles ();
    inline void computeParticleState ();
    inline void collect ();
    inline void initializeSubParticles ();

    // Member declarations

    // Host information
    float* mParticleVertexData;             
    float* mParticleSimulationData;         
    char* mParticleStates;

    // CUDA/OpenGL interoperation information 
    GLuint mParticleVertexDataVbo;
    GLuint mSubParticleVertexDataVbo;
    GLuint mSurfaceParticlesVbo;
    cudaGraphicsResource_t mGraphicsResources[2];
    unsigned int mNumSurfaceParticles;

    //
    // CUDA (device) information
    //
    float* mParticleVertexDataDevPtr;         
    float* mParticleSimulationDataDevPtr;
    char* mParticleStatesDevPtr;
    int* mParticleHashListDevPtr;       
    int* mCellStartListDevPtr;          
    int* mCellEndListDevPtr;         
    int* mIsSurfaceParticleDevPtr;
    unsigned int mThreadsPerBlock;          
    unsigned int mThreadsPerBlockSplit;
    unsigned int mThreadsPerBlockBoundary;
    unsigned int mThreadsPerBlockDefault;
    unsigned int mNumBlocks;                   
    unsigned int mNumBlocksSplit;
    unsigned int mNumBlocksBoundary;
    unsigned int mNumBlocksDefault;

    // Two scale simulation overhead
    int mNumParticlesSplit;
    int mNumParticlesBoundary;
    int mNumParticlesDefault;
    int* mParticleIdsDevPtr;
    int* _isSplitDevPtr;
    int* _isBoundaryDevPtr;
    int* _isDefaultDevPtr;
    int* _splitPrefixSumDevPtr;
    int* _boundaryPrefixSumDevPtr;
    int* _defaultPrefixSumDevPtr;
    int* _subParticleIdListDevPtr;
    int* _defaultParticleIdListDevPtr;
    int* mParticleIdsBoundaryDevPtr;
    int* mParticleIdsSplitDevPtr;
    float* mSubParticleVertexDataDevPtr; 
    float* mSubParticleSimulationDataDevPtr;

    /* Host and device information */
    SimulationParameters _parameters;       /* simulation parameters */

    /* boundary handling information */
    std::string _boundaryMapFileName;
    float* _boundaryMapNodeTableDevPtr;
    unsigned int* _boundaryMapIndexMapDevPtr;

    /* */
    float _leftI;
    float _rightI;

    CudaTimer mTimer;
};

#endif /*include guard of: particle_simulation.h */