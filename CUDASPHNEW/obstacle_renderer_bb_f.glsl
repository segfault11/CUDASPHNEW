#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;
uniform float dX;
uniform float dY;
uniform float dZ;

out vec4 fragColor;

void main() 
{
    fragColor = vec4(0.0f, 0.0f, 1.0f, 1.0f);
	gl_FragDepth = gl_FragCoord.z;
}