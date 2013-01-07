#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;
uniform float particleRadius;

layout(location = 0) in vec3 position;
layout(location = 1) in int particleState;

out VertexData 
{
	float particleRadiusCX;    /* in clip space x direction */
	float particleRadiusCY;    /* in clip space y direction */
	vec4 eye;				   /* eye coordinates of the particle */
    flat int state;
}
vertexData;


void main(void)
{
	/* compute eye coordinates*/
	vertexData.eye = viewMat*vec4(position, 1.0f);
    vertexData.state = particleState & 3;

	/* compute particle radius in clip space in x and y direction
	** NOTE: using screen ratio might save some operations for computing the
	**		 radius in y direction.
	*/
	vertexData.particleRadiusCX = projMat[0][0]*particleRadius/(-vertexData.eye.z);
	vertexData.particleRadiusCY = projMat[1][1]*particleRadius/(-vertexData.eye.z);

	gl_Position = projMat*vertexData.eye/-vertexData.eye.z; /* do perspective division right away */
}