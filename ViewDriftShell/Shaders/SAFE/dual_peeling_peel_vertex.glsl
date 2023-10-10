//--------------------------------------------------------------------------------------
// Order Independent Transparency with Dual Depth Peeling
//
// Author: Louis Bavoil
// Email: sdkfeedback@nvidia.com
//
// Copyright (c) NVIDIA Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

vec3 ShadeVertex() {
	float diffuse = abs( normalize( gl_NormalMatrix * gl_Normal ).z );
	return vec3( gl_Vertex.xy, diffuse );
}


void main( void ) {
	gl_Position = ftransform();
	gl_TexCoord[0].xyz = ShadeVertex();
}
