#include "twoscalestatesub.h"
#include "cgtk\include\opengl.h"
#include <iostream>
using namespace std;

TwoScaleStateSub::TwoScaleStateSub(const ParticleSimulation& sim, 
    unsigned int width, unsigned int height): mSimulation(&sim)
{
    // create glsl object
    mProgram = glCreateProgram();
	cgtkGLAttachShader(mProgram, "twoscalestatesub_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(mProgram, "twoscalestatesub_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLAttachShader(mProgram, "twoscalestatesub_g.glsl", 
        GL_GEOMETRY_SHADER);
	cgtkGLBindFragDataLocation(mProgram, "outFragDepth", 0);
	cgtkGLLinkProgram(mProgram);
	cgtkGLDumpLog(mProgram);

    // initialize glsl object
    mParticleRadius = mSimulation->GetParticleRadius();
    this->setParticleRadius(mParticleRadius);

    mSubParticleRadius = mParticleRadius/2.0f;

    // Init vbo for particle states
    unsigned char* testData = new unsigned char[mSimulation->GetNumParticles()];
    for (unsigned int i = 0; i < mSimulation->GetNumParticles(); i++)
    {
        if (i % 2 == 0)
        {
            testData[i] = 2;
        }
        else
        {
            testData[i] = 3;
        }
    }
    
    //glGenBuffers(1, &_stateVBO);
    //glBindBuffer(GL_ARRAY_BUFFER, _stateVBO);
    //glBufferData(GL_ARRAY_BUFFER, mSimulation->GetNumParticles()*sizeof(unsigned char), 
    //    testData, GL_DYNAMIC_COPY);
    //delete[] testData;

    // initialize vertex array object for base particles
    glGenVertexArrays(1, &mParticleVertexArrayObject);
    glBindVertexArray(mParticleVertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, 
        mSimulation->GetGLParticleVertexBufferObject());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 
        VD_NUM_ELEMENTS*sizeof(float), 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 
        mSimulation->GetGLParticleIndexVertexBufferObject());

    // initialize vertex array object for sub-particles
    glGenVertexArrays(1, &mSubParticleVertexArrayObject);
    glBindVertexArray(mSubParticleVertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, 
        mSimulation->GetGLSubParticleVertexBufferObject());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 
        VD_NUM_ELEMENTS*sizeof(float), 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 
        mSimulation->GetGLSubParticleIndexVertexBufferObject());
}

TwoScaleStateSub::~TwoScaleStateSub()
{
    glDeleteProgram(mProgram);
}

void TwoScaleStateSub::setCamera(float ex, float ey, float ez, float cx, 
	float cy, float cz, float ux, float uy, float uz)
{
		cgtkGLLookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz, mViewMatrix);
		GLint loc = glGetUniformLocation(mProgram, "viewMat");

		if (loc < 0) 
        {
			printf("Could not find location: View Matrix\n");
		}

		glUseProgram(mProgram);
		glUniformMatrix4fv(loc, 1, 0, mViewMatrix);
}

void TwoScaleStateSub::setPerspective(float fovy, float aspect, float n, 
	float f)
{
		cgtkGLPerspective(fovy, aspect, n, f, mProjectionMatrix);
		GLint loc = glGetUniformLocation(mProgram, "projMat");

		if (loc < 0)
        {
			printf("Could not find location: Projection Matrix\n");
		}

		glUseProgram(mProgram);
		glUniformMatrix4fv(loc, 1, 0, mProjectionMatrix);
}

void TwoScaleStateSub::render() const 
{
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(mProgram);

    // draw base particles
    this->setParticleRadius(mParticleRadius);
    glBindVertexArray(mParticleVertexArrayObject);
    glDrawElements(GL_POINTS, mSimulation->GetNumParticlesDefault(), 
        GL_UNSIGNED_INT, 0);

    // draw sub-particles
    this->setParticleRadius(mSubParticleRadius);
    glBindVertexArray(mSubParticleVertexArrayObject);
    glDrawElements(GL_POINTS, mSimulation->GetNumSubParticles(), 
        GL_UNSIGNED_INT, 0);

    glFlush();
	glutSwapBuffers();
	glutPostRedisplay();
}

void TwoScaleStateSub::setParticleRadius (float radius) const
{
    // set particle radius in world space
    GLint loc = glGetUniformLocation(mProgram, "particleRadius");
	glUseProgram(mProgram);
    glUniform1fv(loc, 1, &radius);
}