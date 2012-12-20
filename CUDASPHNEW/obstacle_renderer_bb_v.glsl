#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;
uniform float dX;
uniform float dY;
uniform float dZ;

layout(location = 0) in vec3 position;

void main(void) 
{
	gl_Position = vec4(position, 1.0f);
}