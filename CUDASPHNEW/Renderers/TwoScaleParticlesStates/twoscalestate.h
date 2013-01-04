#ifndef _TWOSCALESTATE_H
#define __TWOSCALESTATE_H

#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include "particle_simulation.h"

/** @class ParticleRenderer
*** @brief Renderer particles as blue spheres.
**/
class TwoScaleState
{
public:
    TwoScaleState(const ParticleSimulation& sim, unsigned int width, 
        unsigned int height);
    ~TwoScaleState();

    void render() const;
    void setCamera(float ex, float ey, float ez, float cx, float cy, float cz,
		float ux, float uy, float uz);
	void setPerspective(float fovy, float aspect, float n, float f);
public:
    static TwoScaleState* example01();

private:
    const ParticleSimulation* _sim;
    GLuint _program;
    GLuint _vertexArrayObject;
    GLfloat _projMat[16];			/* OpenGL style projection matrix */
	GLfloat _viewMat[16];			/* OpenGL style view matrix */
};


#endif /* end of include guard: twoscalestate.h */