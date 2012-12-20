#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;

layout(location = 0) in vec3 position;

void main(void) 
{
	gl_Position = projMat*viewMat*vec4(position, 1.0f);
}