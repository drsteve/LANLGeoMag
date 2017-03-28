///*
// * Vertex shader for Illuminated Lines
// *
// */


uniform vec3 LightPos;
uniform vec3 CameraPos;

attribute vec3 Pprev;
attribute vec3 Pnext;

varying vec3 L, V, T;

void main() {

    vec3 eyeDir;
    vec4 pos;

    //pos    = gl_ModelViewMatrix * gl_Vertex;
    pos  = gl_Vertex;


    T = normalize( Pnext - Pprev );
    V = normalize( CameraPos - pos.xyz );
    L = normalize( LightPos - pos.xyz );

    gl_Position = ftransform();

}

