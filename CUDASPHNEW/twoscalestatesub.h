#ifndef _TWOSCALESTATESUB_H
#define _TWOSCALESTATESUB_H

#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include "particle_simulation.h"

class TwoScaleStateSub
{
public:
    TwoScaleStateSub(const ParticleSimulation& sim, unsigned int width, 
        unsigned int height);
    ~TwoScaleStateSub();

    void render() const;
    void setCamera(float ex, float ey, float ez, float cx, float cy, float cz,
		float ux, float uy, float uz);
	void setPerspective(float fovy, float aspect, float n, float f);
public:
    static TwoScaleStateSub* Example01();

private:
    void setParticleRadius (float radius) const;

    const ParticleSimulation* mSimulation;
    GLuint _stateVBO;
    GLuint mProgram;
    GLuint mParticleVertexArrayObject;
    GLuint mSubParticleVertexArrayObject;
    float mParticleRadius;
    float mSubParticleRadius;
    GLfloat mProjectionMatrix[16];			/* OpenGL style projection matrix */
	GLfloat mViewMatrix[16];			    /* OpenGL style view matrix */
};


#endif /* end of include guard: twoscalestatesub.h */