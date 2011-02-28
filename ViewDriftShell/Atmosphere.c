#include "Atmosphere.h"
#include "Rainbow2.h"

#define Re 6378.135

typedef unsigned char byte;


_aInfo *New_aInfo() {

    _aInfo  *a;
    if ( (a = (_aInfo *) calloc( 1, sizeof( _aInfo ) )) == NULL ) {
        printf("Init_aInfo(): Could not allocate memory for structure\n");
    }
    return( a );

}



/*
   Create an array of sphere vertices. Sphere is centered on origin with radius r, and precision n
   Draw a point for zero radius spheres
*/
void CreateSphereVertices( SphereType *Sphere, int n ) {

    int          i, j, k;
    double       theta1, theta2, Phi, r;
    Lgm_Vector   e, p;

    Sphere->Precision = n;

    r = Sphere->Radius;
    k = 0;

    if (r < 0) r = -r;
    if (n < 0) n = -n;
    if ( (n < 4) || (r <= 0) ) {

        Sphere->nVertex   = 1; // number of vertices needed
        Sphere->Vertex    = (Lgm_Vector *)calloc( Sphere->nVertex, sizeof(Lgm_Vector) );
        Sphere->Color     = (Lgm_Vector *)calloc( Sphere->nVertex, sizeof(Lgm_Vector) );

        Sphere->Vertex[ 0 ].x = 0.0;
        Sphere->Vertex[ 0 ].y = 0.0;
        Sphere->Vertex[ 0 ].z = 0.0;
        Sphere->nVertex = 1;

        return;

    } else {

        Sphere->nVertex   = n*(n+1); // number of vertices needed
        Sphere->Vertex    = (Lgm_Vector *)calloc( Sphere->nVertex, sizeof(Lgm_Vector) );
        Sphere->Color     = (Lgm_Vector *)calloc( Sphere->nVertex, sizeof(Lgm_Vector) );
    }

    for (j=0;j<n/2;j++) {
        // range in latitude from -PI/2 -> PI/2
        theta1 = j * 2.0*M_PI / n - 0.5*M_PI;
        theta2 = (j + 1) * 2.0*M_PI / n - 0.5*M_PI;

//        glBegin(GL_QUAD_STRIP);
        for (i=0;i<=n;i++) {
            Phi = i * 2.0*M_PI / n;


            e.x = cos(theta2) * cos(Phi);
            e.y = cos(theta2) * sin(Phi);
            e.z = sin(theta2);
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;
            //glVertex3f(p.x,p.y,p.z);
            Sphere->Vertex[ k++ ] = p;

            e.x = cos(theta1) * cos(Phi);
            e.y = cos(theta1) * sin(Phi);
            e.z = sin(theta1);
            p.x = r * e.x;
            p.y = r * e.y;
            p.z = r * e.z;
            //glVertex3f(p.x,p.y,p.z);
            Sphere->Vertex[ k++ ] = p;



        }
//        glEnd();
    }

    Sphere->nVertex = k;
    return;

}


/*
   Actually render the vertices as a sphere
*/
void DrawSphereVertices( SphereType *Sphere ) {

    double       Color[3];
    int          i, j, k, n;
    Lgm_Vector   p;

    n = Sphere->Precision;
    k = 0;
    if ( n < 4 ) {
        glDisable( GL_LIGHTING );
        glBegin( GL_POINTS );
        glVertex3f( 0.0, 0.0, 0.0 );
        glEnd( );
        glEnable( GL_LIGHTING );
        return;
    }

    glDisable( GL_LIGHTING );
    for (j=0;j<n/2;j++) {
        glBegin(GL_QUAD_STRIP);
        for (i=0;i<=n;i++) {

            p = Sphere->Vertex[ k++ ];
            if ( Lgm_DotProduct( &aInfo->nCamera, &p ) > -0.3 ) {
                GetScatterColor( p, Color );
                glColor4f( Color[0], Color[1], Color[2], 1.0 );
            } else {
                glColor4f( 0.0, 0.0, 0.0, 1.0 );
            }
            glVertex3f( p.x, p.y, p.z );


            p = Sphere->Vertex[ k++ ];
            if ( Lgm_DotProduct( &aInfo->nCamera, &p ) > -0.3 ) {
                GetScatterColor( p, Color );
                glColor4f( Color[0], Color[1], Color[2], 1.0 );
            } else {
                glColor4f( 0.0, 0.0, 0.0, 1.0 );
            }
            glVertex3f( p.x, p.y, p.z );


        }
        glEnd();
    }
    glEnable( GL_LIGHTING );
    return;

}

void FreeAtmosphere() {

    free( aInfo->DepthBuf );
    free( aInfo->InnerSphere.Vertex );
    free( aInfo->InnerSphere.Color );
    free( aInfo->OuterSphere.Vertex );
    free( aInfo->OuterSphere.Color );
    free( aInfo );

}

void InitAtmosphere() {

    int     j;
    double  u, u2, u3, u4, x, w, w2, w4;

    /*
     *  Set up the Rayleigh scattering parameters.
     */
    aInfo->Kr = 0.0025;     // Rayleigh scattering constant
    aInfo->Kr4PI = aInfo->Kr*4.0*M_PI;
    //aInfo->rScaleHeight = 0.25;
    aInfo->rScaleHeight = 0.45;



    /*
     *  Set up the Mie scattering parameters.
     */
    aInfo->Km = 0.0025;     // Mie scattering constant
    aInfo->Km4PI = aInfo->Km*4.0*M_PI;
    aInfo->mScaleHeight = 0.1;
    aInfo->um = 0.75;       // u parameter needed to compute gm,
                            // the Mie assymmetry factor. See Eq.
                            // (5) of Nishita et al. SIGGRAPH'93
    u = aInfo->um; 
    u2 = u*u; 
    u3 = u2*u; 
    u4 = u2*u2;
    x = 5.0/9.0*u + 125.0/729.0*u3 
        + sqrt( 64.0/27.0 - 325.0/243.0*u2 + 1250.0/2187.0*u4);
    aInfo->gm = 5.0/9.0*u 
        - (4.0/3.0 - 25.0/81.0*u2)*pow(x, -1.0/3.0) 
        + pow(x, 1.0/3.0);
    aInfo->gm2 = aInfo->gm*aInfo->gm;

    /*
     *  Define wavelengths for red, green, and blue channels
     */
    aInfo->SunIntensity  = 15.0;
    aInfo->Sun.x = 1.0; aInfo->Sun.y = 0.0; aInfo->Sun.z = 0.0;
    aInfo->Wavelength[0] = 0.650;   // 650nm Red
    aInfo->Wavelength[1] = 0.570;   // 570nm Grn
    aInfo->Wavelength[2] = 0.475;   // 475nm Blu
    aInfo->nSampleRays   = 6;

    for (j=0; j<3; j++ ) {
        w = aInfo->Wavelength[j];
        w2 = w*w;
        w4 = w2*w2;
        aInfo->Wavelength4[j] = w4;
    }



    /*
     * Define Spheres to represent atmosphere.
     */
    aInfo->OuterSphere.Origin.x = 0.0;
    aInfo->OuterSphere.Origin.y = 0.0;
    aInfo->OuterSphere.Origin.z = 0.0;
    aInfo->OuterSphere.Radius = 1.03;
    aInfo->OuterSphere.Radius2 = aInfo->OuterSphere.Radius*aInfo->OuterSphere.Radius;
    CreateSphereVertices( &aInfo->OuterSphere, 200 );


    aInfo->InnerSphere.Origin.x = 0.0;
    aInfo->InnerSphere.Origin.y = 0.0;
    aInfo->InnerSphere.Origin.z = 0.0;
    aInfo->InnerSphere.Radius = 1.0+10.0/Re;
    aInfo->InnerSphere.Radius2 = aInfo->InnerSphere.Radius*aInfo->InnerSphere.Radius;
    CreateSphereVertices( &aInfo->InnerSphere, 100 );




    aInfo->Scale = 1.0/(aInfo->OuterSphere.Radius - aInfo->InnerSphere.Radius);


    /*
     *  Set parameters for Optical Depth Lookup table
     */
    aInfo->DepthBuf_nSize     = 512;
    aInfo->DepthBuf_nSamples  = 20;
    aInfo->DepthBuf_nChannels = 4;



    /*
     *  Make a lookup table to speed things up
     */
    MakeOpticalDepthBuffer();

}

int iMin(int a, int b) { return( (a<b) ? a : b ); }
int iMax(int a, int b) { return( (a>b) ? a : b ); }

void LookupOpticalDepth( double *p, double x, double y ) {

    double  X, Y, RatioX, RatioY;
    double  *pValue;
    int     nX, nY, i, w, h, nc;

    w = aInfo->DepthBuf_nSize;
    h = aInfo->DepthBuf_nSize;
    nc = aInfo->DepthBuf_nChannels;


    X = x*(w-1); 
    Y = y*(h-1);
    nX = iMin( w-2, iMax(0, (int)X) );
    nY = iMin( h-2, iMax(0, (int)Y) );
    RatioX = X - nX;
    RatioY = Y - nY;

    //pValue = (double *)((unsigned long)aInfo->DepthBuf + nc*sizeof(double) * (w * nY + nX));
    pValue = aInfo->DepthBuf + nc * (w * nY + nX);

    for( i=0; i<nc; i++) {

        p[i] =	pValue[0] * (1.0-RatioX) * (1.0-RatioY) +
                pValue[nc*1] * (RatioX) * (1.0-RatioY) +
                pValue[nc*w] * (1.0-RatioX) * (RatioY) +
                pValue[nc*(w+1)] * (RatioX) * (RatioY);
        pValue++;

    }


}

void MakeOpticalDepthBuffer( ) {

    double      tMinInner, tMinOuter, tMinOuterFwd, tMinInnerFwd;
    double      tMaxInner, tMaxOuter;
    double      rDensityRatio, mDensityRatio;
    double      Cos, Angle, Height, Altitude;
    double      rDepth, mDepth, SampleLength, ScaledLength;
    int         nAngle, nHeight;
    int         nIndex;
    int         i;
    Lgm_Vector  SampleRay;
    RayType     cvRay;
    unsigned char *Image;

    Image = (unsigned char *)calloc( aInfo->DepthBuf_nSize*aInfo->DepthBuf_nSize, sizeof( unsigned char ) );


    aInfo->DepthBuf = (double *)calloc( aInfo->DepthBuf_nSize*aInfo->DepthBuf_nSize*aInfo->DepthBuf_nChannels, sizeof( double) );
    nIndex = 0;
    int iii=0;
// MUST free this somewhere!
	for( nAngle=0; nAngle<aInfo->DepthBuf_nSize; nAngle++) {

		// As the y tex coord goes from 0 to 1, the angle goes from 0 to 180 degrees
		Cos   = 1.0 - (nAngle+nAngle) / (double)aInfo->DepthBuf_nSize;
		Angle = acos( Cos );
        // Ray Direction pointing to the viewpoint
        cvRay.Direction.x = sin( Angle ); cvRay.Direction.y = cos( Angle ); cvRay.Direction.z = 0.0;



		for( nHeight=0; nHeight<aInfo->DepthBuf_nSize; nHeight++ ) {

			// As the x tex coord goes from 0 to 1, the height goes from the bottom of the atmosphere to the top
			Height = 1e-6 + aInfo->InnerSphere.Radius + ((aInfo->OuterSphere.Radius - aInfo->InnerSphere.Radius) * nHeight) / (double)aInfo->DepthBuf_nSize;
            // The position of the camera (origin of ray)
            cvRay.Origin.x = 0.0; cvRay.Origin.y = Height; cvRay.Origin.z = 0.0; 


              
            /*
             *  Do a ray sphere intersect test.  tMinInner is the distance
             *  along the ray from the origin to the first intersection.  If
             *  the ray from vPos heading in the vRay direction intersects the
             *  inner radius (i.e. the planet), then this spot is not visible
             *  from the viewpoint
             */
			if ( !SphereIntersect( &aInfo->InnerSphere, &cvRay, &tMinInner, &tMaxInner, &tMinInnerFwd ) ) {

				rDensityRatio = exp( -(Height - aInfo->InnerSphere.Radius) * aInfo->Scale / aInfo->rScaleHeight );
				mDensityRatio = exp( -(Height - aInfo->InnerSphere.Radius) * aInfo->Scale / aInfo->mScaleHeight );

			} else {

				// Smooth the transition from light to shadow (it is a soft shadow after all)
				rDensityRatio = ((double *)aInfo->DepthBuf)[nIndex - aInfo->DepthBuf_nSize*aInfo->DepthBuf_nChannels] * 0.75;
				mDensityRatio = ((double *)aInfo->DepthBuf)[nIndex+2 - aInfo->DepthBuf_nSize*aInfo->DepthBuf_nChannels] * 0.75;

			}


            /*
             * Determine where the ray intersects the outer radius (the top of
             * the atmosphere) This is the end of our ray for determining the
             * optical depth (cvRay.Origin is the start)
             */
            SphereIntersect( &aInfo->OuterSphere, &cvRay, &tMinOuter, &tMaxOuter, &tMinOuterFwd );


tMinOuterFwd = tMaxOuter;

            /*
             * Next determine the length of each sample, scale the sample ray,
             * and make sure position checks are at the center of a sample ray
             */
			SampleLength = tMinOuterFwd / aInfo->DepthBuf_nSamples;
			ScaledLength = SampleLength * aInfo->Scale;

            SampleRay.x = cvRay.Direction.x * SampleLength;
            SampleRay.y = cvRay.Direction.y * SampleLength;
            SampleRay.z = cvRay.Direction.z * SampleLength;
            cvRay.Origin.x += 0.5*SampleRay.x;
            cvRay.Origin.y += 0.5*SampleRay.y;
            cvRay.Origin.z += 0.5*SampleRay.z;
            /*
             * Iterate through the samples to sum up the optical depth for the
             * distance the ray travels through the atmosphere
             */
			rDepth = 0.0;
			mDepth = 0.0;
			for( i=0; i<aInfo->DepthBuf_nSamples; i++) {

				Height = Lgm_Magnitude( &cvRay.Origin );
				Altitude = (Height - aInfo->InnerSphere.Radius)*aInfo->Scale;
                if (Altitude < 0.0) Altitude = 0.0;
				rDepth += exp( -Altitude / aInfo->rScaleHeight);
				mDepth += exp( -Altitude / aInfo->mScaleHeight);

                cvRay.Origin.x += SampleRay.x;
                cvRay.Origin.y += SampleRay.y;
                cvRay.Origin.z += SampleRay.z;

			}

            /*
             * Multiply the sums by the length the ray traveled
             */
			rDepth *= ScaledLength;
			mDepth *= ScaledLength;




            /*
             * Store the results for Rayleigh to the light source, Rayleigh to
             * the camera, Mie to the light source, and Mie to the camera
             */
            ((double *)aInfo->DepthBuf)[nIndex++] = rDensityRatio;
            ((double *)aInfo->DepthBuf)[nIndex++] = rDepth;
            ((double *)aInfo->DepthBuf)[nIndex++] = mDensityRatio;
            ((double *)aInfo->DepthBuf)[nIndex++] = mDepth;
            *(Image + iii++) = mDensityRatio*1e3;
		}
	}

    FILE *fp_gif = fopen("mDens.gif", "w");
    WriteGIF(fp_gif, (byte *)Image, 0, 512, 512, Rainbow2_Red, Rainbow2_Grn, Rainbow2_Blu, 256, 0, "");
    fclose(fp_gif);



}





// need to input a vertex as an argument.
// the vertex should have position as part of its info
void GetScatterColor( Lgm_Vector v, double Color[3] ){

    int         i, j;
    int         CameraAbove;
    int         CameraInAtmosphere;
    double      tFar, tNear, tMinOuter, tMaxOuter, tMinOuterFwd;
    double      CameraDepth[4], SunDepth[4], SampleDepth[4];
    double      cHeight, cAltitude, vHeight;
    double      SampleLength, ScaledLength;
    double      CameraCosAngle, SunCosAngle, SampleCosAngle;
    double      CosAngle, CosAngle2, rPhase, mPhase;
    double      Height, Altitude;
    double      rDensity, rDepth, mDensity, mDepth;
    double      Atten[4], rSum[3], mSum[3];
    Lgm_Vector  sVec, P;
//    Lgm_Vector  v;      // position of vertex
    Lgm_Vector  c;      // position of camera
    RayType     cvRay;  // ray from camera to vertex



    c =  aInfo->Camera;


    /*
     *  Compute Ray from Camera to Vertex There will be two intersections of
     *  this ray with the sphere. However, since the vertex we are trying to
     *  color is on the sphere, it must be the farthest intersection point
     *  already (i.e. tMaxOuter is just the magnitude of the vector from the
     *  camera to the vertex).
     */
    cvRay.Origin      = c;
    cvRay.Direction.x = v.x - c.x;
    cvRay.Direction.y = v.y - c.y;
    cvRay.Direction.z = v.z - c.z;
    tFar = Lgm_NormalizeVector( &cvRay.Direction );


    /*
     * Perform Ray/Sphere intersection test to get the nearest intersection
     * point.
     */
    SphereIntersect( &aInfo->OuterSphere, &cvRay, &tMinOuter, &tMaxOuter, &tMinOuterFwd );
    tNear = tMinOuterFwd;

    /*
     *  Test to see if we are inside or outside of the atmosphere.
     */
    CameraAbove        = TRUE;
    CameraInAtmosphere = FALSE;
    if ( tNear <= 0.0 ){

        /*
         *  Camera is inside of OuterSphere. I.e., its either inside the
         *  atmosphere (or below the surface of the Earth). 
         */

        // Compute Height and Altitude of Camera.
        CameraInAtmosphere = TRUE;
        tNear = 0.0;
        cHeight   = Lgm_NormalizeVector( &c );                          // Camera Height
        cAltitude = (cHeight - aInfo->InnerSphere.Radius)*aInfo->Scale; // Camera Altitude
        vHeight   = Lgm_Magnitude( &v );                                // Vertex Height
        CameraAbove = ( cHeight >= vHeight );                           // Test if vertex is above camera.
        CameraCosAngle = Lgm_DotProduct( &cvRay.Direction, &c ) / cHeight;
        if ( CameraAbove ) CameraCosAngle *= -1.0;

        /*
         * Use a lookup table to compute the depth from camera to vertex.
         */
        LookupOpticalDepth( CameraDepth, cAltitude, 0.5*(1.0-CameraCosAngle) );


    } else {

        /*
         * We are outside the atmosphere. To avoid integrating over 
         * lots of zeros, shift camera along cvRay right up to where it
         * intersects the outer sphere.
         */
        cvRay.Origin.x += tNear*cvRay.Direction.x;
        cvRay.Origin.y += tNear*cvRay.Direction.y;
        cvRay.Origin.z += tNear*cvRay.Direction.z;
        tFar  -= tNear;
        tNear  = 0.0;
        for (j=0; j<4; j++ ) CameraDepth[j] = 0.0;

    }



 

    /*
     *  If path length of integration is small, bail out -- color is zero
     */
    if ( tFar < MIN_INT_PATH_LEN ) {
        Color[0] = Color[1] = Color[2] = 0.0;
        return;
    }




    /*
     *  Now create loop over a number of sample rays to integrate the scattering
     *  effects into the camera.
     */
    for (j=0; j<3; j++ ){
        rSum[j] = 0.0;
        mSum[j] = 0.0;
    }
    SampleLength = tFar / (double)aInfo->nSampleRays;
    ScaledLength = SampleLength*aInfo->Scale;
    sVec.x = cvRay.Direction.x*SampleLength;
    sVec.y = cvRay.Direction.y*SampleLength;
    sVec.z = cvRay.Direction.z*SampleLength;

    // set position of sample ray
    P.x = cvRay.Origin.x + 0.5*sVec.x; 
    P.y = cvRay.Origin.y + 0.5*sVec.y; 
    P.z = cvRay.Origin.z + 0.5*sVec.z; 
    for ( i=0; i<aInfo->nSampleRays; i++ ){

        Height   = Lgm_Magnitude( &P );
        Altitude = (Height - aInfo->InnerSphere.Radius)*aInfo->Scale;

        /*
         *  Determine optical depth coming from sun to this point.
         */
        SunCosAngle = Lgm_DotProduct( &aInfo->Sun, &P ) / Height;
        LookupOpticalDepth( SunDepth, Altitude, 0.5*(1.0-SunCosAngle) );

        /*
         * Only compute stuff for this sample if SunDepth[0] large enough
         */
        if ( SunDepth[0] > MIN_SUN_DEPTH ) {


            rDensity = ScaledLength * SunDepth[0];
            rDepth   = SunDepth[1];

            mDensity = ScaledLength * SunDepth[2];
            mDepth   = SunDepth[3];


            SampleCosAngle = Lgm_DotProduct( &cvRay.Direction, &P ) / Height;
            if ( CameraAbove ) SampleCosAngle *= -1.0;
            LookupOpticalDepth( SampleDepth, Altitude, 0.5*(1.0-SampleCosAngle) );
            rDepth += (CameraAbove) ? (SampleDepth[1] - CameraDepth[1]) : (CameraDepth[1] - SampleDepth[1]);
            mDepth += (CameraAbove) ? (SampleDepth[3] - CameraDepth[3]) : (CameraDepth[3] - SampleDepth[3]);


            /*
             *  Multiply optical depth by attenuation factor for sample ray
             */
            rDepth *= aInfo->Kr4PI;
            mDepth *= aInfo->Km4PI;


            /*
             *  Compute attenuation factors for sample ray.
             */
            for (j=0; j<3; j++){
                Atten[j] = exp( -rDepth/aInfo->Wavelength4[j] - mDepth );
                rSum[j] += rDensity*Atten[j];
                mSum[j] += mDensity*Atten[j];
            }

        }



        /*
         * Shift location to next sample ray
         */
        P.x += sVec.x;
        P.y += sVec.y;
        P.z += sVec.z;

    }

    /*
     * Calculate the angle between the -cvRay and the Sun
     */
    CosAngle  = -1.0*Lgm_DotProduct( &cvRay.Direction, &aInfo->Sun ); // both ves already normalized
    CosAngle2 = CosAngle*CosAngle;

    rPhase = 0.75*(1.0 + CosAngle2);
    mPhase = 3.0*(1.0 - aInfo->gm2)*(1.0+CosAngle2) / (2.0*(2.0+aInfo->gm2)*pow(1.0+aInfo->gm2-2.0*aInfo->gm*CosAngle, 1.5));

    rPhase *= aInfo->Kr*aInfo->SunIntensity;
    mPhase *= aInfo->Km*aInfo->SunIntensity;

    for ( j=0; j<3; j++ ){
        Color[j] = rSum[j]*rPhase/aInfo->Wavelength4[j] + mSum[j]*mPhase;
        if (Color[j] > 1.0) Color[j] = 1.0;
    }

}

