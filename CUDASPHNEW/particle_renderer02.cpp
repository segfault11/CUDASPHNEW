#include "particle_renderer02.h"
#include "cgtk\include\opengl.h"

ParticleRenderer02::ParticleRenderer02(const ParticleSimulation& sim, 
    unsigned int width, unsigned int height): _sim(&sim)
{
    // create glsl object
    _program = glCreateProgram();
	cgtkGLAttachShader(_program, "particle_renderer02_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(_program, "particle_renderer02_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLAttachShader(_program, "particle_renderer02_g.glsl", 
        GL_GEOMETRY_SHADER);
	cgtkGLBindFragDataLocation(_program, "outFragDepth", 0);
	cgtkGLLinkProgram(_program);
	cgtkGLDumpLog(_program);

    // initialize glsl object

    // set particle radius in world space
    GLint loc = glGetUniformLocation(_program, "particleRadius");
	glUseProgram(_program);
    float radius = _sim->GetParticleRadius();
    glUniform1fv(loc, 1, &radius);

    /*int* test = new int[_sim->GetNumParticles()];
    memset(test, 0, sizeof(int)*_sim->GetNumParticles());*/

    // initialize vertex array object
    glGenVertexArrays(1, &_vertexArrayObject);
    glBindVertexArray(_vertexArrayObject);
    glBindBuffer(GL_ARRAY_BUFFER, _sim->GetGLVertexBufferObject());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 
        VD_NUM_ELEMENTS*sizeof(float), 0);
    glGenBuffers(1, &_isSurfaceParticleVbo);
    glBindBuffer(GL_ARRAY_BUFFER, _isSurfaceParticleVbo);
    glBufferData(GL_ARRAY_BUFFER, _sim->GetNumParticles()*sizeof(int), NULL,
        GL_DYNAMIC_DRAW); 
    glVertexAttribPointer(1, 3, GL_INT, GL_FALSE, sizeof(int), 0);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
}

ParticleRenderer02::~ParticleRenderer02()
{
    glDeleteProgram(_program);
}

void ParticleRenderer02::setCamera(float ex, float ey, float ez, float cx, 
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

void ParticleRenderer02::setPerspective(float fovy, float aspect, float n, 
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

void ParticleRenderer02::render() const 
{
    int* isSurfPartList = 
        ParticleSimulation::CreateIsParticleSurfaceList(_sim);

    glBindBuffer(GL_ARRAY_BUFFER, _isSurfaceParticleVbo);
    glBufferData(GL_ARRAY_BUFFER, _sim->GetNumParticles()*sizeof(int), 
        isSurfPartList, GL_DYNAMIC_DRAW);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindVertexArray(_vertexArrayObject);
    glUseProgram(_program);
    glDrawArrays(GL_POINTS, 0, _sim->GetNumParticles());
    glFlush();
	glutSwapBuffers();
	glutPostRedisplay();

    ParticleSimulation::FreeIsParticleSurfaceList(&isSurfPartList);
}