#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;
uniform float dX;
uniform float dY;
uniform float dZ;

layout (points) in;
layout (line_strip, max_vertices = 24) out;

void main(void) 
{
    
    vec4 p = gl_in[0].gl_Position;
    vec4 pDX = gl_in[0].gl_Position + vec4(dX, 0.0f, 0.0f, 0.0f);
    vec4 pDY = gl_in[0].gl_Position + vec4(0.0f, dY, 0.0f, 0.0f);
    vec4 pDZ = gl_in[0].gl_Position + vec4(0.0f, 0.0f, dZ, 0.0f);
    vec4 pDXDY = gl_in[0].gl_Position + vec4(dX, dY, 0.0f, 0.0f);
    vec4 pDXDZ = gl_in[0].gl_Position + vec4(dX, 0.0f, dZ, 0.0f);
    vec4 pDYDZ = gl_in[0].gl_Position + vec4(0.0f, dY, dZ, 0.0f);
    vec4 pDXDYDZ = gl_in[0].gl_Position + vec4(dX, dY, dZ, 0.0f);

    vec4 pC = projMat*viewMat*p;
    vec4 pDXC = projMat*viewMat*pDX;
    vec4 pDYC = projMat*viewMat*pDY;
    vec4 pDZC = projMat*viewMat*pDZ;
    vec4 pDXDYC = projMat*viewMat*pDXDY;
    vec4 pDXDZC = projMat*viewMat*pDXDZ;
    vec4 pDYDZC = projMat*viewMat*pDYDZ;
    vec4 pDXDYDZC = projMat*viewMat*pDXDYDZ;


    // 1
    gl_Position = pC;
	EmitVertex();

    gl_Position = pDXC;
	EmitVertex();

    EndPrimitive();

    // 2
    gl_Position = pC;
	EmitVertex();

    gl_Position = pDYC;
	EmitVertex();

    EndPrimitive();

    // 3 
    gl_Position = pC;
	EmitVertex();

    gl_Position = pDZC;
	EmitVertex();

    EndPrimitive();

    // 4
    gl_Position = pDXC;
	EmitVertex();

    gl_Position = pDXDYC;
	EmitVertex();

    EndPrimitive();

    // 5
    gl_Position = pDYC;
	EmitVertex();

    gl_Position = pDXDYC;
	EmitVertex();

    EndPrimitive();

    // 6
    gl_Position = pDYC;
	EmitVertex();

    gl_Position = pDYDZC;
	EmitVertex();

    EndPrimitive();

    // 7
    gl_Position = pDZC;
	EmitVertex();

    gl_Position = pDYDZC;
	EmitVertex();

    EndPrimitive();

    // 8
    gl_Position = pDZC;
	EmitVertex();

    gl_Position = pDXDZC;
	EmitVertex();
    
    EndPrimitive();
    
    // 9
    gl_Position = pDYDZC;
	EmitVertex();

    gl_Position = pDXDYDZC;
	EmitVertex();

    EndPrimitive();

    // 10
    gl_Position = pDXDZC;
	EmitVertex();

    gl_Position = pDXDYDZC;
	EmitVertex();

    EndPrimitive();

    // 11
    gl_Position = pDXDYC;
	EmitVertex();

    gl_Position = pDXDYDZC;
	EmitVertex();

    EndPrimitive();

    // 10
    gl_Position = pDXDZC;
	EmitVertex();

    gl_Position = pDXC;
	EmitVertex();

    EndPrimitive();
}