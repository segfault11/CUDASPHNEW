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

/** @class ParticleSimulation
**/
class ParticleSimulation
{
/* public interface */
public:
    /** @brief  Initializes the simulation based on the simulation parameters \
    ***         and the given particle system.
    ***
    *** Allocates memory to store particle information on the device and copies
    *** the initial particle information to it.
    **/
    void init ();

    /** @brief  Makes the Particle Simulation active. 
    ***
    *** Copies the simulation parameters of the particle simulation to the 
    *** device.
    *** Should be used before if another particle simulation was active before.
    **/
    void bind () const;

    /** @brief  Updates the particle positions and attributes in time according \
    ***         to the time step and other parameters set when the constructor  \
    ***         of the particle simulation was invoked.
    **/
    void advance ();

    /** @brief Getter for the OpenGL vertex buffer object, that is used to 
    ***        the particles per vertex information.
    **/
    GLuint getGLVertexBufferObject () const;

    /** @brief Returns the approx. radius of a fluid particle
    **/
    float getParticleRadius () const;
    const char* getParticleState () const;
    unsigned int getNumParticles () const;

    void setNPartThresh (float val);
    void increaseCmDistanceThresh ();
    void decreaseCmDistanceThresh ();

    void saveInfoTable (const std::string& filename);
    void saveParticleInfo (const std::string& filename);


    ~ParticleSimulation ();

    /* static factory methods */

    /** @brief Creates an examples particle simulation.
    ***
    *** In this examples particles for a cube of fluid that is being dropped in 
    *** a rectangular container.
    **/
    static ParticleSimulation* example01 ();
    
    /** @brief Creates a list, that indicates whether a particle belongs to \
    ***        the surface layer or not.
    ***
    *** Creates a list, that indicates whether a particle of simulation belongs 
    *** to the surface layer or not. The list resides in host memory.
    **/
    static int* createIsParticleSurfaceList (const ParticleSimulation* sim);

    /** @brief Frees an isSurfaceParticleList
    **/
    static void freeIsParticleSurfaceList (int** isSurfaceParticleList);

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
    inline void allocateMemoryTwoScale();

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

    /* Member declarations */

    /* Host information */
    float* _particleVertexData;             /* INITIAL vertex data of the
                                            ** particles */
    float* _particleSimulationData;         /* INITIAL simulation data of 
                                            ** the particles */
    char* _particleState;

    // CUDA/OpenGL interoperation information 
    GLuint _particleVbo;
    GLuint _subParticleVbo;
    GLuint _surfaceParticlesVbo;
    cudaGraphicsResource_t _graphicsResource;
    unsigned int _nSurfaceParticles;

    //
    // CUDA (device) information
    //
    float* _particleVertexDataDevPtr;         
    float* _particleSimulationDataDevPtr;
    char* _particleStateDevPtr;
    int* _particleHashListDevPtr;       
    int* _cellStartListDevPtr;          
    int* _cellEndListDevPtr;         
    int* _isSurfaceParticleDevPtr;
    unsigned int _threadsPerBlock;          
    unsigned int _threadsPerBlockSplit;
    unsigned int _threadsPerBlockBoundary;
    unsigned int _threadsPerBlockDefault;
    unsigned int _blocks;                   
    unsigned int _blocksSplit;
    unsigned int _blocksBoundary;
    unsigned int _blocksDefault;


    // Two scale simulation overhead
    int _nSplitParticles;
    int _nBoundaryParticles;
    int _nDefaultParticles;
    int* _particleIdListDevPtr;
    int* _isSplitDevPtr;
    int* _isBoundaryDevPtr;
    int* _isDefaultDevPtr;
    int* _splitPrefixSumDevPtr;
    int* _boundaryPrefixSumDevPtr;
    int* _defaultPrefixSumDevPtr;
    int* _subParticleIdListDevPtr;
    int* _defaultParticleIdListDevPtr;
    int* _boundaryParticleIdListDevPtr;
    int* _splitParticleIdListDevPtr;
    float* _subParticleVertexDataDevPtr; 
    float* _subParticleSimulationDataPtr;

    /* Host and device information */
    SimulationParameters _parameters;       /* simulation parameters */

    /* boundary handling information */
    std::string _boundaryMapFileName;
    float* _boundaryMapNodeTableDevPtr;
    unsigned int* _boundaryMapIndexMapDevPtr;

    /* */
    float _leftI;
    float _rightI;

    CudaTimer _timer;
};

#endif /*include guard of: particle_simulation.h */