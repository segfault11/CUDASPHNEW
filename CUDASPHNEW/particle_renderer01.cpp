#include "particle_renderer01.h"
#include "cgtk\include\opengl.h"

ParticleRenderer::ParticleRenderer(const ParticleSimulation& sim, 
    unsigned int width, unsigned int height): _sim(&sim)
{
    // create glsl object
    _program = glCreateProgram();
	cgtkGLAttachShader(_program, "particle_renderer01_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(_program, "particle_renderer01_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLAttachShader(_program, "particle_renderer01_g.glsl", 
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

    // initialize vertex array object
    glGenVertexArrays(1, &_vertexArrayObject);
    glBindVertexArray(_vertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, _sim->getGLVertexBufferObject());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 
        VD_NUM_ELEMENTS*sizeof(float), 0);
    glEnableVertexAttribArray(0);
}

ParticleRenderer::~ParticleRenderer()
{
    glDeleteProgram(_program);
}

void ParticleRenderer::setCamera(float ex, float ey, float ez, float cx, 
	float cy, float cz, float ux, float uy, float uz)
{
		cgtkGLLookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz, _viewMat);
		GLint loc = glGetUniformLocation(_program, "viewMat");

		if (loc < 0) {
			printf("Could not find location: View Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _viewMat);
}

void ParticleRenderer::setPerspective(float fovy, float aspect, float n, 
	float f)
{
		cgtkGLPerspective(fovy, aspect, n, f, _projMat);
		GLint loc = glGetUniformLocation(_program, "projMat");

		if (loc < 0) {
			printf("Could not find location: Projection Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _projMat);
}

void ParticleRenderer::render() const 
{
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