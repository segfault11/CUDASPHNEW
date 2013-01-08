/*
 * Utilities for Computer Graphics Applications
 * Copyright (C) 2012 Arno Luebke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef OPENGL_H
#define OPENGL_H

#include <GL/GLEW.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*
** OpenGL transformation matrices (View, Perspective). All matrices are 4x4 and
** need to be allocated in Advance. OpenGL matrices are one dimensional arrays
** with 16 (=4*4) elements. Matrices are stored in this array in column-major
** order:
** 				 _               _
** 				| a00 a04 a08 a12 |
** 			M =	| a01 a05 a09 a13 |
** 				| a02 a06 a10 a14 |
** 				|_a03 a07 a11 a15_|
**
** [cgtkGLLookAt] creates the view matrix which transforms from Worldspace to
** Eyespace. ex, ey, ez refer to the position of the camera in world coordinates
** . cx, cy, cz refer to the point the camera looks at. ux, uy, uz refer to the
** vector that points up from the camera's position. mat[16] is a one dim array
** that stores the view matrix. 				 
**
** [cgtkGLFrustum] defines the viewing frustum for the perspective 
** transformation. The input parameters are assumed to be in eye space, and 
** define the viewing frustum. mat[16] stores the projection matrix.
*/
void cgtkGLLookAt(float ex, float ey, float ez, float cx, float cy, float cz,
	float ux, float uy, float uz, GLfloat mat[16]);
void cgtkGLFrustum(float l, float r, float b, float t, float n, float f,
	GLfloat mat[16]);
void cgtkGLPerspective(float fovy, float aspect, float n, float f, 
	GLfloat mat[16]);


void cgtkGLAttachShader(GLuint program, const char* fileName, GLenum type);
void cgtkGLBindAttribLocation(GLuint program, const char* attrName, 
	GLuint index);
void cgtkGLBindFragDataLocation(GLuint program, const char* szColorBufName,
    GLuint index);
void cgtkGLLinkProgram(GLuint program);
void cgtkGLDumpLog(GLuint program);

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: OPENGL_H */
