#include "ViewDriftShell.h"
#include "Vds_DriftShell.h"
typedef struct _MaterialProp {
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess;
} MaterialProp;

static MaterialProp mat_FieldLineType0 = {
// greenish
  {0.0, 0.06, 0.1, 1.0},
  {  0.0/255.0, 105.0/255.0,  27.0/255.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};
static MaterialProp mat_FieldLineType1 = {
// light brownish
  {0.0, 0.06, 0.1, 1.0},
  {152.0/255.0, 153.0/255.0, 122.0/255.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};
static MaterialProp mat_FieldLineType2 = {
// dark brownish
  {0.0, 0.06, 0.1, 1.0},
  { 99.0/255.0*0.6,  99.0/255.0*0.6,  77.0/255.0*0.6, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};
static MaterialProp mat_FieldLineType3 = {
// blueish
  {0.0, 0.06, 0.1, 1.0},
  { 19.0/255.0,  83.0/255.0, 194.0/255.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};
static MaterialProp mat_FieldLineType4 = {
// redish
  {0.0, 0.06, 0.1, 1.0},
  {150.0/255.0,  20.0/255.0,  20.0/255.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};
static MaterialProp mat_FieldLineType5 = {
// yellowish
  {0.0, 0.06, 0.1, 1.0},
  {220.0/255.0, 220.0/255.0,  19.0/255.0, 1.0},
  {0.50196078, 0.50196078, 0.50196078, 1.0},
  0.25
};

/**
 *  \brief
 *      Initialize a Vds_ObjectInfo structure.
 *
 *  \details
 *      Dynamically allocates a Vds_ObjectInfo structure.
 *
 */
Vds_ObjectInfo *Vds_InitObjectInfo( ) {

    Vds_ObjectInfo *o;

    o = (Vds_ObjectInfo *)calloc( 1, sizeof( Vds_ObjectInfo) );
    o->MagEphemAlloced = FALSE;

    return( o );

}


/**
 *  \brief
 *      Free a Vds_ObjectInfo structure.
 *
 *  \details
 *      Frees a dynamically allocated Vds_ObjectInfo structure.
 *
 */
void Vds_FreeObjectInfo( Vds_ObjectInfo *ObjInfo ) {

    if ( ObjInfo->MagEphemAlloced ) {
        free( ObjInfo->MagEphemInfo );
    }
    free( ObjInfo );

    return;

}







// These routines used to be part of gtkglext, but are not included anymore.

/*
 * Quadrics
 */


static GLUquadricObj *lgm_quadObj = NULL;
#define LGM_QUAD_OBJ_INIT() { if(!lgm_quadObj) Lgm_initQuadObj(); }
static void Lgm_initQuadObj( void ) {
  lgm_quadObj = gluNewQuadric();
  if (!lgm_quadObj) g_error("out of memory.");
}

/**
 * Lgm_gl_draw_sphere:
 * @solid: TRUE if the sphere should be solid.
 * @radius: the radius of the sphere.
 * @slices: the number of subdivisions around the Z axis (similar to lines of
 *          longitude).
 * @stacks: the number of subdivisions along the Z axis (similar to lines of
 *          latitude).
 *
 * Renders a sphere centered at the modeling coordinates origin of
 * the specified @radius. The sphere is subdivided around the Z axis into
 * @slices and along the Z axis into @stacks. 
 *
 **/
void Lgm_gl_draw_sphere( gboolean solid, double radius, int slices, int stacks ) {

  LGM_QUAD_OBJ_INIT();

  if (solid) gluQuadricDrawStyle( lgm_quadObj, GLU_FILL );
  else       gluQuadricDrawStyle( lgm_quadObj, GLU_LINE );

  gluQuadricNormals( lgm_quadObj, GLU_SMOOTH );

  /* If we ever changed/used the texture or orientation state
     of quadObj, we'd need to change it to the defaults here
     with gluQuadricTexture and/or gluQuadricOrientation. */
  gluSphere( lgm_quadObj, radius, slices, stacks );

}


/**
 * gdk_gl_draw_cone:
 * @solid: TRUE if the cone should be solid.
 * @base: the radius of the base of the cone.
 * @height: the height of the cone.
 * @slices: the number of subdivisions around the Z axis.
 * @stacks: the number of subdivisions along the Z axis.
 *
 * Renders a cone oriented along the Z axis.
 * The @base of the cone is placed at Z = 0, and the top at Z = @height.
 * The cone is subdivided around the Z axis into @slices, and along
 * the Z axis into @stacks. 
 *
 **/
void Lgm_gl_draw_cone( gboolean solid, double base, double height, int slices, int stacks ) {

  LGM_QUAD_OBJ_INIT();

  if (solid) gluQuadricDrawStyle( lgm_quadObj, GLU_FILL );
  else       gluQuadricDrawStyle( lgm_quadObj, GLU_LINE );

  gluQuadricNormals( lgm_quadObj, GLU_SMOOTH );

  /* If we ever changed/used the texture or orientation state
     of quadObj, we'd need to change it to the defaults here
     with gluQuadricTexture and/or gluQuadricOrientation. */
  gluCylinder( lgm_quadObj, base, 0.0, height, slices, stacks );

}





/**
 *  \brief
 *      Create the field line and drift shell objects (surfaces) from the given
 *      magnetic ephemeris file.
 *
 *  \details
 *      The individual FLs may not all have the same number of points in them.
 *      This makes creating a mesh difficult. An easy solution (adopted here)
 *      is to interpolate tham all to have a prescibed number of points.
 *
 *
 *
 *      \param[in,out]  ObjInfo         Pointer to structure containing object data.
 *
 *
 */
void CreateFieldLinesAndDriftShells( char *Filename, Vds_ObjectInfo *ObjInfo ){

    int          i, tn, ns;


    /*
     *  Read in the MagEphem file
     */
    if ( ObjInfo->MagEphemAlloced ) {
        free( ObjInfo->MagEphemInfo );
        ObjInfo->MagEphemAlloced = FALSE;
    }
    ObjInfo->MagEphemInfo = (Lgm_MagEphemInfo *)calloc( 1, sizeof(Lgm_MagEphemInfo));
    ObjInfo->MagEphemAlloced = TRUE;
    printf( "\t  Reading file: %s\n", Filename );
    ReadMagEphemInfoStruct( Filename, &ObjInfo->nPitchAngles, ObjInfo->MagEphemInfo );
    printf( "\t     Date/Time: %ld/%g\n", ObjInfo->MagEphemInfo->Date, ObjInfo->MagEphemInfo->UTC );
    printf( "\t# Pitch Angles: %d\n", ObjInfo->MagEphemInfo->nAlpha );
if (1==1){
//return;

printf("%d %d %d\n", ObjInfo->nPitchAngles, ObjInfo->MagEphemInfo->nShellPoints[0], ObjInfo->MagEphemInfo->nFieldPnts[0][0]);


    for ( i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) { // Loop over Pitch Angles

        for ( ns=0; ns<ObjInfo->MagEphemInfo->nShellPoints[i]; ns++ ) { // Loop over Shell Field Lines


            /*
             *  Create arrays that include the entire field line from footpoint
             *  to footpoint.
             */
            ObjInfo->nPnts[i][ns] = ObjInfo->MagEphemInfo->nFieldPnts[i][ns];
            for (tn=0; tn<ObjInfo->MagEphemInfo->nFieldPnts[i][ns]; tn++ ){
                ObjInfo->s_gsm[i][ns][tn] = ObjInfo->MagEphemInfo->s_gsm[i][ns][tn];
                ObjInfo->x_gsm[i][ns][tn] = ObjInfo->MagEphemInfo->x_gsm[i][ns][tn];
                ObjInfo->y_gsm[i][ns][tn] = ObjInfo->MagEphemInfo->y_gsm[i][ns][tn];
                ObjInfo->z_gsm[i][ns][tn] = ObjInfo->MagEphemInfo->z_gsm[i][ns][tn];
            }


            /*
             *  Create arrays that include only portion of field between mirror points.
             *  Add the mirror points explicitly later.
             */
///*
int kk;
            for (kk=0, tn=0; tn<ObjInfo->MagEphemInfo->nFieldPnts[i][ns]; tn++ ){
                if ( ( ObjInfo->MagEphemInfo->Bmag[i][ns][tn] < ObjInfo->MagEphemInfo->Bm[i]-1e-6 ) && ( ObjInfo->MagEphemInfo->Bmag[i][ns][tn] > 0.0 ) ) {
                    ObjInfo->s2_gsm[i][ns][kk] = ObjInfo->MagEphemInfo->s_gsm[i][ns][tn];
                    ObjInfo->x2_gsm[i][ns][kk] = ObjInfo->MagEphemInfo->x_gsm[i][ns][tn];
                    ObjInfo->y2_gsm[i][ns][kk] = ObjInfo->MagEphemInfo->y_gsm[i][ns][tn];
                    ObjInfo->z2_gsm[i][ns][kk] = ObjInfo->MagEphemInfo->z_gsm[i][ns][tn];
                    ++kk;
                }
            }
            ObjInfo->nPnts2[i][ns] = kk;
//*/

        }   // Field Line Loop

    }   // Pitch Angle Loop

    MakeFieldLines( 80, ObjInfo ); // FIX -- The 80 should be user setable.
    MakeDriftShellMesh( ObjInfo ); // FIX -- Should be able to control number of points.
}



}


/**
 *  \brief
 *      This routine interpolates the field lines so that each has the same
 *      number of points along their length. We do this in order to be able to
 *      make a regular mesh out of them later.
 *
 *  \details
 *      The individual FLs may not all have the same number of points in them.
 *      This makes creating a mesh difficult. An easy solution (adopted here)
 *      is to interpolate tham all to have a prescibed number of points.
 *
 *
 *
 *      \param[in]      nNewPnts        Number of points desired in the new field lines.
 *      \param[in,out]  ObjInfo         Pointer to structure containing object data.
 *
 *
 */
void MakeFieldLines( int nNewPnts, Vds_ObjectInfo *ObjInfo ) {

    int     i, ns, j;
    double  *xout, *yout, *zout;

    LGM_ARRAY_1D( xout, nNewPnts, double );
    LGM_ARRAY_1D( yout, nNewPnts, double );
    LGM_ARRAY_1D( zout, nNewPnts, double );


    for ( i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++ ) { // Loop over Pitch Angles

        for ( ns=0; ns<ObjInfo->MagEphemInfo->nShellPoints[i]; ns++ ) { // Loop over Filed Lines

            /*
             *  Interp the points to a fixed number so we can make a mesh
             *  out of them.
             */
            if (ObjInfo->nPnts[i][ns] > 4 ){

                InterpFieldLine( ObjInfo->x_gsm[i][ns], ObjInfo->y_gsm[i][ns], ObjInfo->z_gsm[i][ns], ObjInfo->s_gsm[i][ns],
                        ObjInfo->MagEphemInfo->ShellMirror_Ss[i][ns], ObjInfo->MagEphemInfo->ShellMirror_Sn[i][ns],
                        ObjInfo->nPnts[i][ns], xout, yout, zout, nNewPnts );

                for (j=0; j<nNewPnts; j++){
                    ObjInfo->x3_gsm[i][ns][j] = xout[j];
                    ObjInfo->y3_gsm[i][ns][j] = yout[j];
                    ObjInfo->z3_gsm[i][ns][j] = zout[j];
                }

            } else {

                printf("Warning: Line %d in file %s : Insufficient number of field line points: Alpha Index = %d   Field Line Index = %d\n",
                        __LINE__, __FILE__, i, ns );

            }

        }   // Field Line Loop

        if ( ObjInfo->MagEphemInfo->nShellPoints[i] > 0 ) {
            ObjInfo->nFieldPoints[i] = nNewPnts;
        } else {
            ObjInfo->nFieldPoints[i] = 0;
        }

    }   // Pitch Angle Loop

    LGM_ARRAY_1D_FREE( xout );
    LGM_ARRAY_1D_FREE( yout );
    LGM_ARRAY_1D_FREE( zout );


}



/**
 *  \brief
 *      Returns an interpolated version of a field line.
 *
 *  \details
 *      This routine takes a field line defined at n1 points and returns an
 *      interpolated version that has n2 points.
 *
 *
 *
 *      \param[in]          x       Pointer to array containing the x-components of the field line.
 *      \param[in]          y       Pointer to array containing the x-components of the field line.
 *      \param[in]          z       Pointer to array containing the x-components of the field line.
 *      \param[in]         Ss       Distance along the field line of the southern footpoint.
 *      \param[in]         Sn       Distance along the field line of the northern footpoint.
 *      \param[in]         n1       Number of points in the x, y, z arrays.
 *      \param[out]      xout       Pointer to array containing the x-components of the interpolated field line.
 *      \param[out]      yout       Pointer to array containing the x-components of the interpolated field line.
 *      \param[out]      zout       Pointer to array containing the x-components of the interpolated field line.
 *      \param[in]         n1       Number of points in the xout, yout, zout arrays.
 *
 *  \notes
 *      User is responsible for all memory management for the arrays.
 *
 *
 *
 */
void InterpFieldLine( double *x, double *y, double *z, double *s, double Ss, double Sn, int n1, double *xout, double *yout, double *zout, int n2  ){

    gsl_interp_accel    *acc;
    gsl_spline          *spline;
    double  ss, s_inc;
    int     t;

    s_inc = ( Sn-Ss )/((double)(n2-1));


    // Interpolate x positions
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( GSL_INTERP, n1 );
    gsl_spline_init( spline, s, x, n1 );
    //for (ss=s[0]+s_inc, t=0; t<n2; t++ ){
    for (ss=Ss, t=0; t<n2; t++ ){
        xout[t] = gsl_spline_eval( spline, ss, acc );
        ss += s_inc;
    }
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );


    // Interpolate y positions
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( GSL_INTERP, n1 );
    gsl_spline_init( spline, s, y, n1 );
    //for (ss=s[0]+s_inc, t=0; t<n2; t++ ){
    for (ss=Ss, t=0; t<n2; t++ ){
        yout[t] = gsl_spline_eval( spline, ss, acc );
        ss += s_inc;
    }
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );


    // Interpolate z positions
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( GSL_INTERP, n1 );
    gsl_spline_init( spline, s, z, n1 );
    //for (ss=s[0]+s_inc, t=0; t<n2; t++ ){
    for (ss=Ss, t=0; t<n2; t++ ){
        zout[t] = gsl_spline_eval( spline, ss, acc );
        ss += s_inc;
    }
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );

}






void MakeDriftShellMesh( Vds_ObjectInfo *ObjInfo ){

    int                 p, nShellPoints, nAddShellPoints;
    int                 i, j, ns, q, nsp1;
    double              DeltaPhi, PhiInc, Phi, d;
    double              PhiArray[ 2*LGM_LSTARINFO_MAX_FL ], xArray[ 2*LGM_LSTARINFO_MAX_FL ], yArray[ 2*LGM_LSTARINFO_MAX_FL ], zArray[ 2*LGM_LSTARINFO_MAX_FL ];
    double              PhiArray2[ 2*LGM_LSTARINFO_MAX_FL ], xArray2[ 2*LGM_LSTARINFO_MAX_FL ], yArray2[ 2*LGM_LSTARINFO_MAX_FL ], zArray2[ 2*LGM_LSTARINFO_MAX_FL ];
    int                 im1, ip1, k, km1, kp1, MonoTonic;
    Lgm_Vector          u, v, n;
    gsl_interp_accel    *acc_x, *acc_y, *acc_z;
    gsl_spline          *spline_x, *spline_y, *spline_z;


    /*
     *  We typically have more than enough points all the driftshell field
     *  lines.  However, we usually dont have very many field lines. So the
     *  drift shell surface can look more faceted than we would like. To solve
     *  this, we want to add more points in-between the drift shell points.
     *
     *  Here, we assume that we have already interped the FL points to give a
     *  constant number of points along each FL.  If we take one index number
     *  from each FL, the points should describe a closed curve. Its these
     *  curves that we can interpolate on to get a more smooth overall surface.
     */
//why is this hard-coded?
//nShellPoints = 24;
//why is this hard-coded?
nAddShellPoints = 10; // number of points to add in between each given Shell point
    for (p=0; p<ObjInfo->MagEphemInfo->nAlpha; p++) { // Loop over pitch angle

        nShellPoints = ObjInfo->MagEphemInfo->nShellPoints[p];

        for (i=0; i<ObjInfo->nFieldPoints[p]; i++ ){

            // create interpolation arrays
            for (ns=0; ns<nShellPoints; ns++ ){

                /*
                 * Compute Phi -- the azimuthal angle of the points.
                 */
                PhiArray[ns] = DegPerRad*atan2( ObjInfo->y3_gsm[p][ns][i],  ObjInfo->x3_gsm[p][ns][i] );
                xArray[ns] = ObjInfo->x3_gsm[p][ns][i];
                yArray[ns] = ObjInfo->y3_gsm[p][ns][i];
                zArray[ns] = ObjInfo->z3_gsm[p][ns][i];

            }

            /*
             * Make sure the Phi Array does not have values that are the same.
             */
            for (MonoTonic = TRUE, ns=1; ns<nShellPoints; ns++ ){
                if ( fabs( PhiArray[ns] - PhiArray[ns-1]) < 1e-3 ) {
                    MonoTonic = FALSE;
                }
            }

            // KLUDGE
            if ( !MonoTonic ) {
                for ( ns=0; ns<nShellPoints; ns++ ){
                    PhiArray[ns] = ns*360.0/(double)nShellPoints;
                }
            }


            // keep first phi val as is.
            PhiArray2[0] = PhiArray[0];
            xArray2[0]   = xArray[0];
            yArray2[0]   = yArray[0];
            zArray2[0]   = zArray[0];
            d = 0.0;
            for (ns=1; ns<nShellPoints; ns++ ){
                if ( PhiArray[ns] < PhiArray[ns-1] ) d += 360.0;
                PhiArray2[ns] = PhiArray[ns]+d;
                xArray2[ns] = xArray[ns];
                yArray2[ns] = yArray[ns];
                zArray2[ns] = zArray[ns];
            }
            /*
             * Finally Add in the first point again for periodic spline
             */
            PhiArray2[nShellPoints] = PhiArray[0]+d;
            xArray2[nShellPoints] = xArray[0];
            yArray2[nShellPoints] = yArray[0];
            zArray2[nShellPoints] = zArray[0];
            if ( PhiArray2[nShellPoints] < PhiArray[nShellPoints-1] ) PhiArray2[nShellPoints] += 360.0;






            /*
             *  Add 'nAddShellPoints' points between each of the nominal points. This makes the surface smoother.
             */
            acc_x    = gsl_interp_accel_alloc( ); spline_x = gsl_spline_alloc( gsl_interp_cspline_periodic, nShellPoints+1 ); gsl_spline_init( spline_x, PhiArray2, xArray2, nShellPoints+1 );
            acc_y    = gsl_interp_accel_alloc( ); spline_y = gsl_spline_alloc( gsl_interp_cspline_periodic, nShellPoints+1 ); gsl_spline_init( spline_y, PhiArray2, yArray2, nShellPoints+1 );
            acc_z    = gsl_interp_accel_alloc( ); spline_z = gsl_spline_alloc( gsl_interp_cspline_periodic, nShellPoints+1 ); gsl_spline_init( spline_z, PhiArray2, zArray2, nShellPoints+1 );
            for (j=0, ns=0; ns<nShellPoints; ns++ ){

                nsp1 = ns+1; // last value is at index nShellPoints (so ns+1 is allowed)
                DeltaPhi = PhiArray2[nsp1] - PhiArray2[ns];
                PhiInc = DeltaPhi/((double)(nAddShellPoints+1));

                // add original point
                ObjInfo->x4_gsm[p][j][i] = xArray2[ns];
                ObjInfo->y4_gsm[p][j][i] = yArray2[ns];
                ObjInfo->z4_gsm[p][j][i] = zArray2[ns];
                ++j; // count point added so far

                // add new points
                for (q=0; q<nAddShellPoints; q++){
                    Phi = PhiArray2[ns] + PhiInc*(double)(q+1);
                    ObjInfo->x4_gsm[p][j][i] = gsl_spline_eval( spline_x, Phi, acc_x );
                    ObjInfo->y4_gsm[p][j][i] = gsl_spline_eval( spline_y, Phi, acc_y );
                    ObjInfo->z4_gsm[p][j][i] = gsl_spline_eval( spline_z, Phi, acc_z );
                    ++j; // count point added so far
                }

                if (ns==nShellPoints-1){
                    // add first point again to close
                    ObjInfo->x4_gsm[p][j][i] = xArray2[0];
                    ObjInfo->y4_gsm[p][j][i] = yArray2[0];
                    ObjInfo->z4_gsm[p][j][i] = zArray2[0];
                    ++j; // count point added so far
                }
            }
            gsl_spline_free( spline_x ); gsl_interp_accel_free( acc_x );
            gsl_spline_free( spline_y ); gsl_interp_accel_free( acc_y );
            gsl_spline_free( spline_z ); gsl_interp_accel_free( acc_z );
            ObjInfo->nShellPoints4 = j;

        }


        /*
         *  Compute Normals
         */
        for (i=0;i<ObjInfo->nFieldPoints[p];i++) {
            for (k=0;k<ObjInfo->nShellPoints4;k++) {

                if (( i>0 ) && ( i<(ObjInfo->nFieldPoints[p]-1) )) {
                    /*
                     * Interior FL points
                     */
                    im1 = i-1; ip1 = i+1;
                } else if ( i == 0 ) {
                    /*
                     * Start FL point
                     */
                    im1 = i; ip1 = i+1;
                } else {
                    /*
                     * End FL point
                     */
                    im1 = i-1; ip1 = i;
                }
                u.x = ObjInfo->x4_gsm[p][k][ip1] - ObjInfo->x4_gsm[p][k][im1];
                u.y = ObjInfo->y4_gsm[p][k][ip1] - ObjInfo->y4_gsm[p][k][im1];
                u.z = ObjInfo->z4_gsm[p][k][ip1] - ObjInfo->z4_gsm[p][k][im1];

                km1 = k-1; kp1 = k+1;
                if (km1 < 0) km1 += ObjInfo->nShellPoints4;
                if (kp1 >= ObjInfo->nShellPoints4) kp1 -= ObjInfo->nShellPoints4;
                v.x = ObjInfo->x4_gsm[p][kp1][i] - ObjInfo->x4_gsm[p][km1][i];
                v.y = ObjInfo->y4_gsm[p][kp1][i] - ObjInfo->y4_gsm[p][km1][i];
                v.z = ObjInfo->z4_gsm[p][kp1][i] - ObjInfo->z4_gsm[p][km1][i];

                Lgm_CrossProduct( &u, &v, &n );
                Lgm_NormalizeVector( &n );
                ObjInfo->nx4_gsm[p][k][i] = -n.x;
                ObjInfo->ny4_gsm[p][k][i] = -n.y;
                ObjInfo->nz4_gsm[p][k][i] = -n.z;

            }
        }

    } // p loop


}



void GenerateDriftShellLists( Vds_ObjectInfo *ObjInfo ){

    int i, p, k;


    /*
     *  Create List for the Drift Shell Surfaces
     */
    ObjInfo->DriftShellList4 = glGenLists( ObjInfo->MagEphemInfo->nAlpha );

    // Drift Shell Surfaces
    for (p=0; p<ObjInfo->MagEphemInfo->nAlpha; p++){

        glNewList( ObjInfo->DriftShellList4 + p, GL_COMPILE );
            for (i=0;i<ObjInfo->nFieldPoints[p]-1;i++) {
                glBegin(GL_QUAD_STRIP);
                for (k=0;k<ObjInfo->nShellPoints4;k++) {
                //for (k=0;k<ObjInfo->nShellPoints4*3/4;k++) {
                    glNormal3f( ObjInfo->nx4_gsm[p][k][i], ObjInfo->ny4_gsm[p][k][i], ObjInfo->nz4_gsm[p][k][i] );
                    glVertex3f( ObjInfo->x4_gsm[p][k][i], ObjInfo->y4_gsm[p][k][i], ObjInfo->z4_gsm[p][k][i] );
                    glNormal3f( ObjInfo->nx4_gsm[p][k][i+1], ObjInfo->ny4_gsm[p][k][i+1], ObjInfo->nz4_gsm[p][k][i+1] );
                    glVertex3f( ObjInfo->x4_gsm[p][k][i+1], ObjInfo->y4_gsm[p][k][i+1], ObjInfo->z4_gsm[p][k][i+1] );
                }
                glEnd();
            }
        glEndList( );

    }

}
void ReGenerateDriftShellLists( Vds_ObjectInfo *ObjInfo ){
    glDeleteLists( ObjInfo->DriftShellList4, ObjInfo->MagEphemInfo->nAlpha );
    GenerateDriftShellLists( ObjInfo );
}


void GenerateFieldLineLists( Vds_ObjectInfo *ObjInfo ){

    int i, ns;

if (1==1){
// dont think we get these anymore...
    /*
     *  Create List for the Full Field Lines
     */
    ObjInfo->DriftShellList2 = glGenLists( ObjInfo->MagEphemInfo->nAlpha );

    // Field Lines and Foot Points
    for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++){

        if ( ObjInfo->MagEphemInfo->Lstar[i] > 1.0 ) {

            glNewList( ObjInfo->DriftShellList2 + i, GL_COMPILE );

            //glMaterialfv( GL_FRONT, GL_DIFFUSE, colors[i] );
            for (ns=0; ns<ObjInfo->MagEphemInfo->nShellPoints[i]; ns++){

                // North Foot Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();

                // South Foot Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();


                // Field Lines -- What is the diff between these two lines here?
MakeTube( ObjInfo->x_gsm[i][ns], ObjInfo->y_gsm[i][ns], ObjInfo->z_gsm[i][ns], ObjInfo->nPnts[i][ns], 12, 0.0375/4.0 );
//int ppp;
//for ( ppp=0; ppp<ObjInfo->nPnts[i][ns]; ++ppp){
//printf("ObjInfo->x_gsm[i][ns][ppp], ObjInfo->z_gsm[i][ns][ppp] = %g %g\n", ObjInfo->x_gsm[i][ns][ppp], ObjInfo->z_gsm[i][ns][ppp]);
//}
//exit(0);
//                MakeTube( ObjInfo->x2_gsm[i][ns], ObjInfo->y2_gsm[i][ns], ObjInfo->z2_gsm[i][ns], ObjInfo->nPnts2[i][ns], 12, 0.0375/2.0 );
//                MakeTube( ObjInfo->x3_gsm[i][ns], ObjInfo->y3_gsm[i][ns], ObjInfo->z3_gsm[i][ns], ObjInfo->nFieldPoints[i], 12, 0.0375/2.0 );


                // North Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();

                // South Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();

            }
            glEndList( );

        }
    }

}




    /*
     *  Create List for the Partial Field Lines
     *  (I think this is all we get now out of the MagEphem files... Correct?)
     */
    ObjInfo->DriftShellList3 = glGenLists( ObjInfo->MagEphemInfo->nAlpha );

    // Field Lines and Foot Points
    for (i=0; i<ObjInfo->MagEphemInfo->nAlpha; i++){
        glNewList( ObjInfo->DriftShellList3 + i, GL_COMPILE );

        for (ns=0; ns<ObjInfo->MagEphemInfo->nShellPoints[i]; ns++){
        //for (ns=0; ns<ObjInfo->MagEphemInfo->nShellPoints[i]*3/4; ns++){

                // North Foot Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();


                // South Foot Points
                glPushMatrix(); glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();



                // Field Lines
//printf("ObjInfo->nFieldPoints[i] = %d\n", ObjInfo->nFieldPoints[i]);
//exit(0);
                MakeTube( ObjInfo->x3_gsm[i][ns], ObjInfo->y3_gsm[i][ns], ObjInfo->z3_gsm[i][ns], ObjInfo->nFieldPoints[i], 12, 0.0375/2.0 );



                // North Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();


                // South Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].z ); 
                Lgm_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();

            }



        glEndList( );
    }

}


void ReGenerateFieldLineLists( Vds_ObjectInfo *ObjInfo ){
    glDeleteLists( ObjInfo->DriftShellList2, ObjInfo->MagEphemInfo->nAlpha );
    glDeleteLists( ObjInfo->DriftShellList3, ObjInfo->MagEphemInfo->nAlpha );
    GenerateFieldLineLists( ObjInfo );
}


/*
 *   Kludge to do arbitrary FLs
 */
void GenerateMiscFieldLineLists( Vds_ObjectInfo *ObjInfo ){

    double  x, y, z;
    int     i, j, ns, Gap, Flag, Type;
    char    Line[128];
    char    Line2[128];
    FILE    *fp;

    fp = fopen( "/home/mgh//git/LanlGeoMag/Examples/Trace/FieldLines.txt", "r" );
LGM_ARRAY_1D( ObjInfo->nPnts5,   2000, int );
LGM_ARRAY_1D( ObjInfo->FL_Type5, 2000, int );
LGM_ARRAY_2D( ObjInfo->x5_gsm,   2000, 2000, double );
LGM_ARRAY_2D( ObjInfo->y5_gsm,   2000, 2000, double );
LGM_ARRAY_2D( ObjInfo->z5_gsm,   2000, 2000, double );

    j = -1; i = 0; Flag = 0;
    while ( fgets( Line, 80, fp ) != EOF ) {

        if ( Line[0] == '\0' ) {
            break;
        } else if ( Line[0] == 'T' ) {

            // new FL
            ++j;
            sscanf( Line, "Type: %d\n", &Type );
            ObjInfo->FL_Type5[j] = Type;
            //printf("=========\n");
            //printf("ObjInfo->FL_Type5[%d] = %d\n", j, ObjInfo->FL_Type5[j] );

            i = 0;

        } else {

            sscanf( Line, "%lf %lf %lf %d\n", &x, &y, &z, &Gap );
            //printf("j = %d   i=%d   x, y, z, Gap = %g %g %g %d\n", j, i, x, y, z, Gap);
            ObjInfo->x5_gsm[j][i] = x;
            ObjInfo->y5_gsm[j][i] = y;
            ObjInfo->z5_gsm[j][i] = z;
            if (i<1999) ++i;
            ObjInfo->nPnts5[j] = i;
        }



        Line[0]='\0';
    }
//exit(0);
    fclose( fp );
    ObjInfo->nFLs = j;


    /*
     *  Create List for the Full Field Lines
     */
    ObjInfo->MiscFieldLines = glGenLists( 1 );
    glNewList( ObjInfo->MiscFieldLines, GL_COMPILE );
    for (j=0; j<ObjInfo->nFLs; j++){
        switch ( ObjInfo->FL_Type5[j] ){
        case 0:
            // OPEN
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType3.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType3.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType3.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType3.shininess * 128.0);
            break;
        case 1:
            // CLOSED DAY
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType0.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType0.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType0.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType0.shininess * 128.0);
            break;
        case 2:
            // CLOSED NIGHT y>20
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType1.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType1.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType1.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType1.shininess * 128.0);
            break;
        case 3:
            // CLOSED NIGHT 
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType2.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType2.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType2.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType2.shininess * 128.0);
            break;
        case 4:
            // LOBE
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType3.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType3.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType3.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType3.shininess * 128.0);
            break;
        case 5:
            // NORTH FUNNY
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType4.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType4.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType4.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType4.shininess * 128.0);
            break;
        case 6:
            // SOUTH FUNNY
            glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_FieldLineType5.ambient);
            glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_FieldLineType5.diffuse);
            glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_FieldLineType5.specular);
            glMaterialf(  GL_FRONT, GL_SHININESS, mat_FieldLineType5.shininess * 128.0);
            break;
        }

        MakeTube( ObjInfo->x5_gsm[j], ObjInfo->y5_gsm[j], ObjInfo->z5_gsm[j], ObjInfo->nPnts5[j], 6, 3*0.1/7 );
        //MakeTube( ObjInfo->x5_gsm[j], ObjInfo->y5_gsm[j], ObjInfo->z5_gsm[j], ObjInfo->nPnts5[j], 12, 0.3 );
    }
//exit(0);
    glEndList( );



}

void ReGenerateMiscFieldLineLists( Vds_ObjectInfo *ObjInfo ){
    glDeleteLists( ObjInfo->MiscFieldLines, 1 );
LGM_ARRAY_1D_FREE( ObjInfo->nPnts5 );
LGM_ARRAY_1D_FREE( ObjInfo->FL_Type5 );
LGM_ARRAY_2D_FREE( ObjInfo->x5_gsm );
LGM_ARRAY_2D_FREE( ObjInfo->y5_gsm );
LGM_ARRAY_2D_FREE( ObjInfo->z5_gsm );
    GenerateMiscFieldLineLists( ObjInfo );
}



