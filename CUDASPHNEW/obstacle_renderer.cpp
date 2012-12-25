#include "obstacle_renderer.h"
#include "cgtk\include\opengl.h"
#include <stdio.h>

ObstacleRenderer::ObstacleRenderer() 
{
    // create glsl object
    _program = glCreateProgram();
	cgtkGLAttachShader(_program, "obstacle_renderer_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(_program, "obstacle_renderer_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLBindFragDataLocation(_program, "fragColor", 0);
	cgtkGLLinkProgram(_program);
	cgtkGLDumpLog(_program);

    // create glsl object
    _program2 = glCreateProgram();
	cgtkGLAttachShader(_program2, "obstacle_renderer_bb_v.glsl", 
        GL_VERTEX_SHADER);
	cgtkGLAttachShader(_program2, "obstacle_renderer_bb_f.glsl", 
        GL_FRAGMENT_SHADER);
	cgtkGLAttachShader(_program2, "obstacle_renderer_bb_g.glsl", 
        GL_GEOMETRY_SHADER);
	cgtkGLBindFragDataLocation(_program2, "fragColor", 0);
	cgtkGLLinkProgram(_program2);
	cgtkGLDumpLog(_program2);




    // set up opengl buffers
    glGenVertexArrays(1, &_vao);
    glBindVertexArray(_vao);
    glGenBuffers(1, &_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glGenBuffers(1, &_indexVbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _indexVbo);

    // set up opengl buffers
    glGenVertexArrays(1, &_vao2);
    glBindVertexArray(_vao2);
    glGenBuffers(1, &_vbo2);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo2);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);


}

ObstacleRenderer::~ObstacleRenderer() 
{
    glDeleteProgram(_program);
    glDeleteVertexArrays(1, &_vao);
    glDeleteBuffers(1, &_vbo);
}

void ObstacleRenderer::setCamera(float ex, float ey, float ez, float cx, 
	float cy, float cz, float ux, float uy, float uz)
{
		cgtkGLLookAt(ex, ey, ez, cx, cy, cz, ux, uy, uz, _viewMat);
		GLint loc = glGetUniformLocation(_program, "viewMat");
        GLint loc2 = glGetUniformLocation(_program2, "viewMat");

		if (loc < 0) {
			printf("Could not find location: View Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _viewMat);
        glUseProgram(_program2);
		glUniformMatrix4fv(loc2, 1, 0, _viewMat);
}

void ObstacleRenderer::setPerspective(float fovy, float aspect, float n, 
	float f)
{
		cgtkGLPerspective(fovy, aspect, n, f, _projMat);
		GLint loc = glGetUniformLocation(_program, "projMat");
        GLint loc2 = glGetUniformLocation(_program2, "projMat");

		if (loc < 0) {
			printf("Could not find location: Projection Matrix\n");
		}

		glUseProgram(_program);
		glUniformMatrix4fv(loc, 1, 0, _projMat);

        glUseProgram(_program2);
		glUniformMatrix4fv(loc2, 1, 0, _projMat);
}

void ObstacleRenderer::setObstacle(const TriangleMesh& obs)
{
    _obstacle = &obs;

    glBindVertexArray(_vao);

    // update indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _indexVbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
        sizeof(unsigned int)*3*_obstacle->_nFaces, _obstacle->_faceList,
        GL_STATIC_DRAW);

    // update vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*_obstacle->_nVertices*3, 
        _obstacle->_vertexList, GL_STATIC_DRAW);





    Rectangle3f bb = obs.getBoundingBox();

    bb.dump();

    float or[3]; // origin of the bounding box
    or[0] = bb.getV1().getX();
    or[1] = bb.getV1().getY();
    or[2] = bb.getV1().getZ();
    float dx = bb.getV2().getX() - or[0]; 
    float dy = bb.getV2().getY() - or[1]; 
    float dz = bb.getV2().getZ() - or[2]; 

    Vector3f(dx,dy,dz).dump();


    glBindVertexArray(_vao2);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo2);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3, or, GL_STATIC_DRAW);

    GLint loc = glGetUniformLocation(_program2, "dX");
    if (loc < 0) {
	    printf("Could not find location: dX\n");
    }
    glUniform1f(loc, dx);

    loc = glGetUniformLocation(_program2, "dY");
    if (loc < 0) {
	    printf("Could not find location: dY\n");
    }
    glUniform1f(loc, dy);

    loc = glGetUniformLocation(_program2, "dZ");
    if (loc < 0) {
	    printf("Could not find location: dZ\n");
    }
    glUniform1f(loc, dz);


}

void ObstacleRenderer::draw() const
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnableVertexAttribArray(0);
    

    glUseProgram(_program2);
	glBindVertexArray(_vao2);
    glDrawArrays(GL_POINTS, 0, 1);



    glBindVertexArray(_vao);
    glUseProgram(_program);
    glDrawElements(GL_TRIANGLES, 3*_obstacle->_nFaces, 
        GL_UNSIGNED_INT, 0); 

    glFlush();
    glutSwapBuffers();
    glutPostRedisplay();
}