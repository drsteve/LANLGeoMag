///*
// * Vertex shader for Illuminated Lines
// *
// */


uniform vec3 LightPos;
uniform vec3 CameraPos;

attribute vec3 Tang;

varying vec3 L, V, T;

void main() {

    vec4 pos;
    pos  = gl_Vertex;

    T = Tang;
    V = normalize( CameraPos - pos.xyz );
    L = normalize( LightPos - pos.xyz );

    gl_Position = ftransform();

}

