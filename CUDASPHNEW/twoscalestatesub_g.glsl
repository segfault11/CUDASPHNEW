#version 330 core
precision highp float;

uniform float particleRadius; /* in world space */
uniform uint colorIndex;
uniform vec3 particleColor;

layout (points) in;

in VertexData 
{
	float particleRadiusCX;
	float particleRadiusCY;
	vec4 eye;
}
vertexData[1];

layout (triangle_strip, max_vertices = 4) out;

out GeometryData 
{
	vec4 eye;
	vec2 relCoord;
} 
geometryData;


void main (void) 
{
    float dx = vertexData[0].particleRadiusCX;
	float dy = vertexData[0].particleRadiusCY;

	gl_Position = gl_in[0].gl_Position + vec4(+dx, +dy, 0.0f, 0.0f);
	geometryData.eye = vertexData[0].eye;
	geometryData.relCoord = vec2(1.0f, 1.0f);
	EmitVertex();
	
	gl_Position = gl_in[0].gl_Position + vec4(-dx, +dy, 0.0f, 0.0f);
	geometryData.eye = vertexData[0].eye;
	geometryData.relCoord = vec2(-1.0f, 1.0f);
	EmitVertex();

	gl_Position = gl_in[0].gl_Position + vec4(+dx, -dy, 0.0f, 0.0f);
	geometryData.eye = vertexData[0].eye;
	geometryData.relCoord = vec2(1.0f, -1.0f);
	EmitVertex();

	gl_Position = gl_in[0].gl_Position + vec4(-dx, -dy, 0.0f, 0.0f);
	geometryData.eye = vertexData[0].eye;
	geometryData.relCoord = vec2(-1.0f, -1.0f);
 	EmitVertex();

    EndPrimitive();
}