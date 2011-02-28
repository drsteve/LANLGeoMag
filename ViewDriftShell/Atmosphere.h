#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Lgm/Lgm_Vec.h>
#include "Objects.h"
#include <GL/gl.h>
#include <GL/glu.h>




#define MIN_INT_PATH_LEN    1e-6
#define MIN_SUN_DEPTH       1e-6

typedef struct _aInfo {

    /*
     *  Rayleigh scattering parameters
     */
    double  Kr;
    double  Kr4PI;
    double  rScaleHeight;


    /*
     *  Mie scattering parameters
     */
    double  Km;
    double  Km4PI;
    double  um;
    double  gm;
    double  gm2;
    double  mScaleHeight;

    /*
     *  Wavelengths for Red,Grn,Blu
     */
    double  Wavelength[3];
    double  Wavelength4[3];
    double  SunIntensity;


    /*
     * Positions of the Sun and the Camera
     */
    Lgm_Vector  Sun;
    Lgm_Vector  Camera;
    Lgm_Vector  nCamera; // normalized
    Lgm_Vector  Up;
    Lgm_Vector  Right;
    double      Rcam; // magnitude of Camera.


    /*
     * Spheres defining atmosphere.
     */
    SphereType   OuterSphere;
    SphereType   InnerSphere;
    double       Scale;


    /*
     * Number of sample rays to use in the scattering calcs
     */
    int nSampleRays;


    /*
     * Depth Buffer params
     */
    int     DepthBuf_nSize;
    int     DepthBuf_nSamples;
    int     DepthBuf_nChannels;
    double  *DepthBuf;





} _aInfo;


_aInfo  *New_aInfo();
void    InitAtmosphere();
void    FreeAtmosphere();
int     iMin(int a, int b);
int     iMax(int a, int b);
void    LookupOpticalDepth( double *p, double x, double y );
void    MakeOpticalDepthBuffer( );
void    GetScatterColor( Lgm_Vector v, double Color[3] );
void    CreateSphereVertices( SphereType *Sphere, int n );
void    DrawSphereVertices( SphereType *Sphere );


/*
 * Global defs
 */
#ifdef MAIN
/*
 *  Make aInfo Global. We could get around this,
 *  but it adds complexity.
 */
_aInfo   *aInfo;
#else
extern _aInfo   *aInfo;
#endif




#endif

