#ifndef _OBSTACLE_RENDERER_H
#define _OBSTACLE_RENDERER_H

#include <Windows.h>
#include <gl\glew.h>
#include <gl\glut.h>
#include "obstacle.h"

class ObstacleRenderer
{
public:
    ObstacleRenderer();
    ~ObstacleRenderer();

    void setCamera(float ex, float ey, float ez, float cx, 
	    float cy, float cz, float ux, float uy, float uz);
    void setPerspective(float fovy, float aspect, float n, 
	    float f);
    void setObstacle(const Obstacle& obs);
    void draw() const;

private:
    GLuint _program;
    GLuint _program2;
    GLuint _vao;
    GLuint _vbo;
    GLuint _vao2;
    GLuint _vbo2;
    GLuint _indexVbo;
    GLfloat _projMat[16];			/* OpenGL style projection matrix */
	GLfloat _viewMat[16];			/* OpenGL style view matrix */
    const Obstacle* _obstacle;
};

#endif /* end of include guard: obstacle_renderer.h */ 