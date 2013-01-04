#include "twoscalestate.h"
#include "cgtk\include\opengl.h"
#include <iostream>
using namespace std;

TwoScaleState::TwoScaleState(const ParticleSimulation& sim, 
    unsigned int width, unsigned int height): _sim(&sim)
{
    // create glsl object
    _program = glCreateProgram();
	cgtkGLAttachShader(_program, "twoscalestate_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(_program, "twoscalestate_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLAttachShader(_program, "twoscalestate_g.glsl", 
        GL_GEOMETRY_SHADER);
	cgtkGLBindFragDataLocation(_program, "outFragDepth", 0);
	cgtkGLLinkProgram(_program);
	cgtkGLDumpLog(_program);

    // initialize glsl object

    // set particle radius in world space
    GLint loc = glGetUniformLocation(_program, "particleRadius");
	glUseProgram(_program);
    float radius = _sim->getParticleRadius();
    glUniform1fv(loc, 1, &radius);

    // init vbo for particle states
    unsigned char* testData = new unsigned char[_sim->getNumParticles()];
    for (unsigned int i = 0; i < _sim->getNumParticles(); i++)
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
    
    glGenBuffers(1, &_stateVBO);
    glBindBuffer(GL_ARRAY_BUFFER, _stateVBO);
    glBufferData(GL_ARRAY_BUFFER, _sim->getNumParticles()*sizeof(unsigned char), 
        testData, GL_DYNAMIC_COPY);
    delete[] testData;

    // initialize vertex array object
    glGenVertexArrays(1, &_vertexArrayObject);
    glBindVertexArray(_vertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, _sim->getGLVertexBufferObject());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 
        VD_NUM_ELEMENTS*sizeof(float), 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, _stateVBO);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_BYTE, 0, 0);
    glEnableVertexAttribArray(1);
    cout << hex << glGetError() << endl;
    
}

TwoScaleState::~TwoScaleState()
{
    glDeleteProgram(_program);
}

void TwoScaleState::setCamera(float ex, float ey, float ez, float cx, 
	float cy, float cz, float ux, float uy, float uz)
{
		cgtkGLLookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz, _viewMat);
		GLint loc = glGetUniformLocation(_program, "viewMat");

		if (loc < 0) 
        {
			printf("Could not find location: View Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _viewMat);
}

void TwoScaleState::setPerspective(float fovy, float aspect, float n, 
	float f)
{
		cgtkGLPerspective(fovy, aspect, n, f, _projMat);
		GLint loc = glGetUniformLocation(_program, "projMat");

		if (loc < 0)
        {
			printf("Could not find location: Projection Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _projMat);
}

void TwoScaleState::render() const 
{
    glBindBuffer(GL_ARRAY_BUFFER, _stateVBO);
    glBufferData(GL_ARRAY_BUFFER, _sim->getNumParticles()*sizeof(unsigned char), 
        _sim->getParticleState(), GL_DYNAMIC_COPY);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindVertexArray(_vertexArrayObject);
    glUseProgram(_program);
    glDrawArrays(GL_POINTS, 0, _sim->getNumParticles());
    glFlush();
	glutSwapBuffers();
	glutPostRedisplay();
}