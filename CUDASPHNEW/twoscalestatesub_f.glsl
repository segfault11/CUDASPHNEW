#version 330 core
precision highp float;

uniform mat4 projMat;
uniform mat4 viewMat;
uniform float particleRadius; /* in world space */
uniform vec3 particleColor;

in GeometryData 
{
	vec4 eye;
	vec2 relCoord;	
} 
geometryData;

out vec4 outFragDepth;

void main(void) 
{
	/* compute the normal of the sphere at the current fragment location 
	** discard, if fragment isnt part of the sphere
	*/
	vec4 n = vec4(geometryData.relCoord, 0.0f, 1.0f);

	float norm = n.x*n.x + n.y*n.y;

	if (norm > 1.0f) 
    {
		discard;
	}

	n.z = sqrt(1 - norm);

    vec4 colour;
    colour = vec4(0.0f, 0.0f, 1.0f, 1.0f);

	/* compute eye-space position of the fragment using n
	** NOTE: for this to work eye.z should be normalized by eye.w.
	**       at the moment this is not necessary as there are no
	**		 object transformations (rotation, translation etc).
	*/
	float zEye = geometryData.eye.z + particleRadius*n.z;

	/* compute the depth values of this fragment in NDC */
	float z = -projMat[2][2] - projMat[2][3]/zEye;

	vec4 l = vec4(0.0, 0.0, 1.0, 0.0f);
	float diffuse = max(0.0, dot(n, l));
    float atten = 0.9 + 0.1*diffuse;

	//float scale = (geometryData.energy)/(320.0f);
	//outFragDepth = diffuse*scale*vec4(1.0f, 0.0f, 0.0f, 1.0f);
	outFragDepth = atten*vec4(particleColor, 1.0f);

	gl_FragDepth = z;
}