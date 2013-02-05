#include "ViewDriftShell.h"
#include "Vds_DriftShell.h"


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

    int     i, tn, ns;


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
/*
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
*/

        }   // Field Line Loop

    }   // Pitch Angle Loop

    MakeFieldLines( 80, ObjInfo ); // FIX -- The 80 should be user setable.
    MakeDriftShellMesh( ObjInfo ); // FIX -- Should be able to control number of points.

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

        ObjInfo->nFieldPoints[i] = nNewPnts;

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
    double              PhiArray[50], xArray[50], yArray[50], zArray[50];
    double              PhiArray2[50], xArray2[50], yArray2[50], zArray2[50];
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
nShellPoints = 24;
//why is this hard-coded?
nAddShellPoints = 10; // number of points to add in between each given Shell point
    for (p=0; p<ObjInfo->MagEphemInfo->nAlpha; p++) { // Loop over pitch angle
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

if (0==1){
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
                gdk_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();

                // South Foot Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();


                // Field Lines -- What is the diff between these two lines here?
                MakeTube( ObjInfo->x2_gsm[i][ns], ObjInfo->y2_gsm[i][ns], ObjInfo->z2_gsm[i][ns], ObjInfo->nPnts2[i][ns], 12, 0.0375/2.0 );
//                MakeTube( ObjInfo->x3_gsm[i][ns], ObjInfo->y3_gsm[i][ns], ObjInfo->z3_gsm[i][ns], ObjInfo->nFieldPoints[i], 12, 0.0375/2.0 );


                // North Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();

                // South Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
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

                // North Foot Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Pn[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();


                // South Foot Points
                glPushMatrix(); glTranslatef(ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellSphericalFootprint_Ps[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.01, 30, 30 ); 
                glPopMatrix();



                // Field Lines
                MakeTube( ObjInfo->x3_gsm[i][ns], ObjInfo->y3_gsm[i][ns], ObjInfo->z3_gsm[i][ns], ObjInfo->nFieldPoints[i], 12, 0.0375/2.0 );



                // North Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Pn[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
                glPopMatrix();


                // South Mirror Points
                glPushMatrix(); 
                glTranslatef(ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].x, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].y, 
                             ObjInfo->MagEphemInfo->ShellMirror_Ps[i][ns].z ); 
                gdk_gl_draw_sphere( TRUE, 0.10/4.0, 30, 30 ); 
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



