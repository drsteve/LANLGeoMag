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

uniform vec3 LightPos;
uniform vec3 CameraPos;

#define MAX_DEPTH 1.0

varying vec3 N;
varying vec3 v;    


vec4 ShadeFragment( void ) {  

    //vec3 L = normalize(gl_LightSource[0].position.xyz - v);   
    //vec3 E = normalize(vec3(10,0,0)-v); // we are in Eye Coordinates, so EyePos is (0,0,0)  
    //vec3 R = normalize(-reflect(L,N));  

    vec3 L = normalize( LightPos - v );   
    vec3 E = normalize( CameraPos - v ); 
    vec3 R = normalize( -reflect( L,normalize(N) ) );  
 
    //calculate Ambient Term:  
    vec4 Iamb = gl_FrontLightProduct[0].ambient;    

    //calculate Diffuse Term:  
    //vec4 Idiff = gl_FrontLightProduct[0].diffuse * max(dot(N,L), 0.0);
    vec4 Idiff = gl_FrontLightProduct[0].diffuse * abs(dot(N,L));
    //vec4 Idiff = gl_FrontLightProduct[0].diffuse * dot(N,L);
    Idiff = clamp(Idiff, 0.0, 1.0);     
   
    // calculate Specular Term:
    vec4 Ispec = gl_FrontLightProduct[0].specular * pow(dot(R,E),gl_FrontMaterial.shininess);
    Ispec = clamp(Ispec, 0.0, 1.0); 

    // write Total Color:  
    vec4 Itotal = gl_FrontLightModelProduct.sceneColor + Iamb + Idiff + Ispec;
    //Itotal.a = gl_FrontLightProduct[0].diffuse.a;


    return Itotal;
}
          

void main(void) {

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
color.a = 1.00;
	gl_FragData[0].xy = vec2(-MAX_DEPTH);
	
	if (fragDepth == nearestDepth) {
		gl_FragData[1].xyz += color.rgb * color.a * alphaMultiplier;
		gl_FragData[1].w = 1.0 - alphaMultiplier * (1.0 - color.a);
	} else {
		gl_FragData[2] += color;
	}
}
