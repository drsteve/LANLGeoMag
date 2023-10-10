//--------------------------------------------------------------------------------------
// Order Independent Transparency with Dual Depth Peeling
//
// Author: Louis Bavoil
// Email: sdkfeedback@nvidia.com
//
// Copyright (c) NVIDIA Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

varying vec3 N;
varying vec3 v;

void main(void) {     

    //v = vec3(gl_ModelViewMatrix * gl_Vertex);       
    //N = normalize(gl_NormalMatrix * gl_Normal);
    v = gl_Vertex;       
    N = gl_Normal;

    //gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;  
    gl_Position = ftransform();

}

