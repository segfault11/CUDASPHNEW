#ifndef _PARTICLE_RENDERER_H
#define _PARTICLE_RENDERER_H

#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include "particle_simulation.h"

/** @class ParticleRenderer
*** @brief Renderer particles as blue spheres.
**/
class ParticleRenderer
{
public:
    ParticleRenderer(const ParticleSimulation& sim, unsigned int width, 
        unsigned int height);
    ~ParticleRenderer();
    void render() const;
    void setCamera(float ex, float ey, float ez, float cx, float cy, float cz,
		float ux, float uy, float uz);
	void setPerspective(float fovy, float aspect, float n, float f);
public:
    static ParticleRenderer* example01();

private:
    const ParticleSimulation* _sim;
    GLuint _program;
    GLuint _vertexArrayObject;
    GLfloat _projMat[16];			/* OpenGL style projection matrix */
	GLfloat _viewMat[16];			/* OpenGL style view matrix */
};


#endif /* end of include guard: particle_renderer.h */