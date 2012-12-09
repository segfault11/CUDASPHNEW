#include "../include/opengl.h"
#include "../include/math.h"
#include "../include/error.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

char* read_file(const char* filename);

void cgtkGLLookAt(float ex, float ey, float ez, float cx, float cy, float cz,
	float ux, float uy, float uz, GLfloat mat[16]) {
	CGTKVec3f n, up, u, v;

	/* compute base vectors of the camera space, i.e. u, v, n */
	n[0] = cx - ex;
	n[1] = cy - ey;
	n[2] = cz - ez;
	cgtkVec3fNormalize(n);

	up[0] = ux;
	up[1] = uy;
	up[2] = uz;

	cgtkVec3fCrossProduct(n, up, u);
	cgtkVec3fNormalize(u);

	cgtkVec3fCrossProduct(u, n, v);

	/* fill transformation matrix */
	mat[0] = u[0];
	mat[1] = v[0];
	mat[2] = -n[0];
	mat[3] = 0.0f;
	mat[4] = u[1];
	mat[5] = v[1];
	mat[6] = -n[1];
	mat[7] = 0.0f;
	mat[8] = u[2];
	mat[9] = v[2];
	mat[10] = -n[2];
	mat[11] = 0.0f;
	mat[12] = -u[0]*ex - u[1]*ey - u[2]*ez;
	mat[13] = -v[0]*ex - v[1]*ey - v[2]*ez;
	mat[14] = n[0]*ex + n[1]*ey + n[2]*ez;
	mat[15] = 1.0f;
}

void cgtkGLFrustum(float l, float r, float b, float t, float n, float f,
	GLfloat mat[16])
{
	mat[0] = 2.0f*n/(r - l);
	mat[1] = 0.0f;
	mat[2] = 0.0f;
	mat[3] = 0.0f;
	
	mat[4] = 0.0f;
	mat[5] = 2.0f*n/(t - b);
	mat[6] = 0.0f;
	mat[7] = 0.0f;

	mat[8] = (r + l)/(r - l);
	mat[9] = (t + b)/(t - b);
	mat[10] = (f + n)/(n - f);
	mat[11] = -1.0f;

	mat[12] = 0.0f;
	mat[13] = 0.0f;
	mat[14] = -2.0f*f*n/(f - n);
	mat[15] = 0.0f;
}

void cgtkGLPerspective(float fovy, float aspect, float n, float f, 
	GLfloat mat[16]) {
		float t = tanf((fovy/2.0)*(3.141593f/180));
		float h = n*t;
		float w = h*aspect;
		cgtkGLFrustum(-w, w, -h, h, n, f, mat);
}

void cgtkGLAttachShader(GLuint program, const char* filename, GLenum type) {
	GLuint shader;
	GLint hasCompiled;
	char* fileContents;	

	if (program == 0) {
	    CGTK_DUMP_ERROR("Invalid program handle.");
		return;
	}

	/* read contents froms file */
	fileContents = read_file(filename);

	if (fileContents == NULL) {
	    CGTK_DUMP_ERROR("File not found.");
		return;
	}

	/* create shader */
	shader = glCreateShader(type);

	if (shader == 0) {
		CGTK_DUMP_ERROR("Could not create shader");
		free(fileContents);
		return;
	}
	
	/* attach source code to shader object */
	glShaderSource(shader, 1, (const GLchar**) &fileContents, NULL);
	
	/* compile shader and check for errors */
	glCompileShader(shader);
    free(fileContents);
	glGetShaderiv(shader, GL_COMPILE_STATUS, &hasCompiled);
    
    if (!hasCompiled) {
		CGTK_DUMP_ERROR("Could not compile shader");
		glDeleteShader(shader);
		return;
	}

	/* attach shader to the program */
	glAttachShader(program, shader);
}

/*
** Associates a generic vertex attribute index with a named attribute variable.
*/
void cgtkGLBindAttribLocation(GLuint program, const char* attrName, 
	GLuint index) {
	glBindAttribLocation(program, index, attrName);
}

void cgtkGLBindFragDataLocation(GLuint program, const char* szColorBufName,
    GLuint index) {
	glBindFragDataLocation(program, index, szColorBufName);
}

void cgtkGLLinkProgram(GLuint program) {
	GLint hasLinked;
	
	glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &hasLinked);

    if (!hasLinked) 
    {
		CGTK_DUMP_ERROR("Could not link program.");
		return;
    }
}

void cgtkGLDumpLog(GLuint program) {
	GLint logLength;
	char* log;
	
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
	log = (char*) malloc((logLength + 1) * sizeof(char));
	glGetProgramInfoLog(program, logLength, NULL, log);
	log[logLength] = '\0';
	
	printf("%s\n", log);
	
	free(log); 
}

char* read_file(const char* filename) {
    FILE* pFile;
    long fileSize = 0;
    char* pContents;
    char ch;

    pFile = fopen(filename, "r");

    if (pFile == NULL) 
    {
        return NULL;
    }

    while (1) {
        ch = fgetc(pFile);
        if (ch == EOF)
            break;
        ++fileSize;
    }
	rewind (pFile);
    pContents = (char*) malloc((fileSize + 1)*sizeof(char));

    if (pContents == NULL) 
    {
        fclose(pFile);
        return NULL;
    }

    fread(pContents, sizeof(char), fileSize, pFile);
    pContents[fileSize] = '\0';
    fclose(pFile);

    return pContents;
}
