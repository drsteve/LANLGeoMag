//--------------------------------------------------------------------------------------
// Order Independent Transparency with Dual Depth Peeling
//
// Author: Louis Bavoil
// Email: sdkfeedback@nvidia.com
//
// Copyright (c) NVIDIA Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

#extension ARB_draw_buffers : require

uniform samplerRECT DepthBlenderTex;
uniform samplerRECT FrontBlenderTex;

#define MAX_DEPTH 1.0




uniform float Alpha;

#define COLOR_FREQ 30.0
#define ALPHA_FREQ 30.0

vec4 ShadeFragment() {

	float xWorldPos = gl_TexCoord[0].x;
	float yWorldPos = gl_TexCoord[0].y;
	float diffuse = gl_TexCoord[0].z;

	vec4 color;
	float i = floor(xWorldPos * COLOR_FREQ);
	float j = floor(yWorldPos * ALPHA_FREQ);
	color.rgb = (fmod(i, 2.0) == 0) ? vec3(.4,.85,.0) : vec3(1.0);
	//color.a = (fmod(j, 2.0) == 0) ? Alpha : 0.2;
    Alpha = 1.0; // mgh
	color.a = Alpha;


	color.rgb *= diffuse;
	return color;
}


void main( void ) {

	// window-space depth interpolated linearly in screen space
	float fragDepth = gl_FragCoord.z;

	vec2 depthBlender = textureRect(DepthBlenderTex, gl_FragCoord.xy).xy;
	vec4 forwardTemp = textureRect(FrontBlenderTex, gl_FragCoord.xy);
	
	// Depths and 1.0-alphaMult always increase
	// so we can use pass-through by default with MAX blending
	gl_FragData[0].xy = depthBlender;
	
	// Front colors always increase (DST += SRC*ALPHA_MULT)
	// so we can use pass-through by default with MAX blending
	gl_FragData[1] = forwardTemp;
	
	// Because over blending makes color increase or decrease,
	// we cannot pass-through by default.
	// Each pass, only one fragment writes a color greater than 0
	gl_FragData[2] = vec4(0.0);

	float nearestDepth = -depthBlender.x;
	float farthestDepth = depthBlender.y;
	float alphaMultiplier = 1.0 - forwardTemp.w;

	if (fragDepth < nearestDepth || fragDepth > farthestDepth) {
		// Skip this depth in the peeling algorithm
		gl_FragData[0].xy = vec2(-MAX_DEPTH);
		return;
	}
	
	if (fragDepth > nearestDepth && fragDepth < farthestDepth) {
		// This fragment needs to be peeled again
		gl_FragData[0].xy = vec2(-fragDepth, fragDepth);
		return;
	}
	
	// If we made it here, this fragment is on the peeled layer from last pass
	// therefore, we need to shade it, and make sure it is not peeled any farther
	vec4 color = ShadeFragment();
	gl_FragData[0].xy = vec2(-MAX_DEPTH);
	
	if (fragDepth == nearestDepth) {
		gl_FragData[1].xyz += color.rgb * color.a * alphaMultiplier;
		gl_FragData[1].w = 1.0 - alphaMultiplier * (1.0 - color.a);
	} else {
		gl_FragData[2] += color;
	}

}
