///*
// * Fragment shader for Illuminated Lines
// *
// *  The fragment shader needs to receive L, V, T and compute B, N, H as follows
// *      First re-normalize L, V, T
// *
// *      B = T cross V
// *      N = B cross T
// *      H = (V+L)/|V+L|
// *
// *  From these we need to compute;
// *
// *      LT       = L dot T
// *      CosAlpha = (L dot N)/sqrt(1-(L dot T)^2) 
// *      CosBeta  = (H dot N)/sqrt(1-(H dot T)^2) 
// *
// *  (Note , that quantities coming in from the vertex shader will be
// *  interpolated to the position of the current fragment. This means that the
// *  vectors L, V, T will not be properly normalized. So they must be normalized
// *  each time in the frag shader -- first step above.)
// *
// *  Frag shader also needs access to the two 2D texture images that we need to
// *  precompute in the main code.  These texture images are; Fd( CosAlpha, LT )
// *  and Fs( CosAlpha, CosBeta ). These are just intensity images -- e.g. store
// *  with only one componeber, e.g. like GL_RED image format. Once the Frag
// *  shader has all of thes values, the final color out is given by;
// *
// *      fd = Fd( CosAlpha, LT )
// *      fs = Fs( CosAlpha, CosBeta )
// *      Cout = Cin*( ka + kd*fd + ks*sqrt( 1-HT*HT )*fs )
// */
//
uniform sampler2D FdTexture; // a GL_RED image (only one component)
uniform sampler2D FsTexture; // a GL_RED image (only one component)

//uniform vec4 LineColor; // Base color of line
uniform float ka;       // ambient factor
uniform float kd;       // diffuse factor
uniform float ks;       // specular factor
uniform float n;        // specular exponent

varying vec3 L, V, T;
varying vec4 C;

void main() {

    float   fd, fs, LT, HT, CosAlpha, CosBeta, SqrtOmHT2;
    vec3    B, H, N, Ln, Vn, Tn, ColorOut;
    vec2    TexCoords;

    Ln = normalize( L );
    Vn = normalize( V );
    Tn = normalize( T );

    H = normalize( Vn + Ln );

    B = normalize( cross( Tn, Vn ) );
    N = normalize( cross( B, Tn ) );

    LT  = dot( Ln, Tn );
    HT  = dot( H, Tn );
    SqrtOmHT2 = sqrt(1.0 - HT*HT);

    CosAlpha = dot( Ln, N ) / sqrt( 1.0 - LT*LT );
    CosBeta  = dot( H,  N ) / SqrtOmHT2;

    float tCosAlpha = (CosAlpha+1.0)/2.0;
    float tCosBeta  = (CosBeta+1.0)/2.0;
    float tLT       = (LT+1.0)/2.0;
    

    /*
     * Retrieve values of fd and fs from the textures.
     * texture() returns a vec4. The four compoments (s,t,p,q) can be swizzled.
     * We only need to select the first component since our textures should
     * have only one component.
     */
    TexCoords = vec2( tCosAlpha, tLT );
    fd = texture2D( FdTexture, TexCoords ).s; 

    TexCoords = vec2( tCosAlpha, tCosBeta );
    fs = texture2D( FsTexture, TexCoords ).s;


    ColorOut = C.rgb*( ka + kd*fd + ks*pow( SqrtOmHT2, n )*fs );
    gl_FragColor = vec4( ColorOut, C.a );

}
