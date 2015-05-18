#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Lgm/Lgm_Vec.h"
#include "Lgm/Lgm_HDF5.h"
#include "Lgm/Lgm_JPLeph.h"

#ifndef LGM_INDEX_DATA_DIR
#warning "hard-coding LGM_INDEX_DATA_DIR because it was not in config.h"
#define LGM_INDEX_DATA_DIR    /usr/local/share/LanlGeoMag/Data
#endif


Lgm_JPLephemInfo *Lgm_InitJPLephemInfo( int DEnum, int getBodies, int verbosity ) {

    Lgm_JPLephemInfo  *jpl = (Lgm_JPLephemInfo *) calloc (1, sizeof(Lgm_JPLephemInfo));
    Lgm_InitJPLephDefaults( DEnum, getBodies, verbosity, jpl );
    return jpl;

}

void Lgm_InitJPLephDefaults( int DEnum, int getBodies, int verbosity, Lgm_JPLephemInfo *jpl ) {

    //getBodies is (logical OR) sum of variable numbers
    //Sun=1; EarthMoon=2; InnerPlanets=4; OuterPlanets=8; LibrationNutation=16;
    jpl->DEnum = DEnum;
    jpl->verbosity = verbosity;

    if (LGM_DE_SUN & getBodies) { jpl->getSun=TRUE; }
    if (LGM_DE_EARTHMOON & getBodies) { jpl->getEarth=TRUE; }
    if (LGM_DE_INNERPLANETS & getBodies) { jpl->getInnerPlanets=TRUE; }
    if (LGM_DE_OUTERPLANETS & getBodies) { jpl->getOuterPlanets=TRUE; }
    if (LGM_DE_LIBR_NUT & getBodies) { jpl->getLibrationNutation=TRUE; }

}

void Lgm_FreeJPLephemInfo( Lgm_JPLephemInfo *jpl ) {

    if ( jpl->SunAlloced ) { 
        LGM_ARRAY_3D_FREE( jpl->sun );
    }

    if ( jpl->EarthMoonAlloced ) {
        LGM_ARRAY_3D_FREE( jpl->earthmoon );
        LGM_ARRAY_3D_FREE( jpl->moon_wrt_earth );
    }

    if ( jpl->InnerPlanetsAlloced ) {
        LGM_ARRAY_3D_FREE( jpl->mercury );
        LGM_ARRAY_3D_FREE( jpl->venus );
        LGM_ARRAY_3D_FREE( jpl->mars );
    }

    if ( jpl->OuterPlanetsAlloced ) {
        LGM_ARRAY_3D_FREE( jpl->jupiter );
        LGM_ARRAY_3D_FREE( jpl->saturn );
        LGM_ARRAY_3D_FREE( jpl->uranus );
        LGM_ARRAY_3D_FREE( jpl->neptune );
        LGM_ARRAY_3D_FREE( jpl->pluto );
    }

    if ( jpl->LibrNutAlloced ) {
        LGM_ARRAY_3D_FREE( jpl->libration );
        LGM_ARRAY_3D_FREE( jpl->nutation );
    }

    free( jpl );

}


void Lgm_ReadJPLephem( Lgm_JPLephemInfo *jpl ) {

    double      *JDparams;
    char        *Path, JPLephemPath[2048];
    char        eph_num[14], JPLephemFile[2048];
    int         StatError, InFileExists;
    struct stat StatBuf;
    herr_t      status;
    hid_t       file;
    hsize_t     dims[4];

    /*
     * read from HDF5 file
     */
    Path = getenv( "JPL_EPHEM_PATH" );
    if ( Path == NULL ) {

        strcpy( JPLephemPath, LGM_INDEX_DATA_DIR );
        strcat( JPLephemPath, "/DE_FILES" );

    } else {
        /*
         * Test for existence
         */
        struct stat sts;
        if ( ( stat( Path, &sts ) ) == -1 ) {
            strcpy( JPLephemPath, LGM_INDEX_DATA_DIR );
            strcat( JPLephemPath, "/DE_FILES" );
            printf("Environment variable JPL_EPHEM_PATH points to a non-existent directory: %s. Setting JPLephemPath to: %s \n", Path, JPLephemPath );
        } else {
            strcpy( JPLephemPath, Path );
        }

    }

    // jpl structure has member DEnum that should be cast to a string
    sprintf( eph_num, "/jpl_de%d.h5", jpl->DEnum );
    strcpy( JPLephemFile, JPLephemPath );
    strcat( JPLephemFile, eph_num );

    InFileExists = ( (StatError = stat( JPLephemFile, &StatBuf )) != -1 ) ? TRUE : FALSE;

    if ( ( InFileExists ) && ( H5Fis_hdf5( JPLephemFile ) )) {

        // Read HDF5 file here
        file = H5Fopen( JPLephemFile, H5F_ACC_RDONLY, H5P_DEFAULT );
        if (jpl->verbosity > 1) {
            printf("Loading JPL definitive ephemeris from %s\n", JPLephemFile);
        }
        // get DEXXX Chebyshev polynomial coefficients for requested objects read into variables
        if (jpl->getSun) {
            jpl->sun = Get_DoubleDataset_3D( file, "/Sun", dims);
            jpl->sun_nvals = dims[0]; jpl->sun_naxes = dims[1]; jpl->sun_ncoeffs = dims[2];
            jpl->SunAlloced = TRUE;
        }
        if (jpl->getEarth) { //gets Earth-Moon barycenter (heliocentric) and Moon (geocentric)
            jpl->earthmoon = Get_DoubleDataset_3D( file, "/EarthMoon", dims);
            jpl->earthmoon_nvals = dims[0];
            jpl->earthmoon_naxes = dims[1];
            jpl->earthmoon_ncoeffs = dims[2];
            jpl->moon_wrt_earth = Get_DoubleDataset_3D( file, "/Moon", dims);
            jpl->moon_wrt_earth_nvals = dims[0];
            jpl->moon_wrt_earth_naxes = dims[1];
            jpl->moon_wrt_earth_ncoeffs = dims[2];
            jpl->EarthMoonAlloced = TRUE;
        }
        if (jpl->getInnerPlanets) { 
            jpl->mercury = Get_DoubleDataset_3D( file, "/Mercury", dims);
            jpl->mercury_nvals = dims[0]; jpl->mercury_naxes = dims[1]; jpl->mercury_ncoeffs = dims[2];
            jpl->venus = Get_DoubleDataset_3D( file, "/Venus", dims);
            jpl->venus_nvals = dims[0]; jpl->venus_naxes = dims[1]; jpl->venus_ncoeffs = dims[2];
            jpl->mars = Get_DoubleDataset_3D( file, "/Mars", dims);
            jpl->mars_nvals = dims[0]; jpl->mars_naxes = dims[1]; jpl->mars_ncoeffs = dims[2];
            jpl->InnerPlanetsAlloced = TRUE;
        }
        if (jpl->getOuterPlanets) { 
            jpl->jupiter = Get_DoubleDataset_3D( file, "/Jupiter", dims);
            jpl->jupiter_nvals = dims[0]; jpl->jupiter_naxes = dims[1]; jpl->jupiter_ncoeffs = dims[2];
            jpl->saturn = Get_DoubleDataset_3D( file, "/Saturn", dims);
            jpl->saturn_nvals = dims[0]; jpl->saturn_naxes = dims[1]; jpl->saturn_ncoeffs = dims[2];
            jpl->uranus = Get_DoubleDataset_3D( file, "/Uranus", dims);
            jpl->uranus_nvals = dims[0]; jpl->uranus_naxes = dims[1]; jpl->uranus_ncoeffs = dims[2];
            jpl->neptune = Get_DoubleDataset_3D( file, "/Neptune", dims);
            jpl->neptune_nvals = dims[0]; jpl->neptune_naxes = dims[1]; jpl->neptune_ncoeffs = dims[2];
            jpl->pluto = Get_DoubleDataset_3D( file, "/Pluto", dims);
            jpl->pluto_nvals = dims[0]; jpl->pluto_naxes = dims[1]; jpl->pluto_ncoeffs = dims[2];
            jpl->OuterPlanetsAlloced = TRUE;
        }
        if (jpl->getLibrationNutation) { 
            jpl->libration = Get_DoubleDataset_3D( file, "/Librations", dims);
            jpl->libration_nvals = dims[0]; jpl->libration_naxes = dims[1]; jpl->libration_ncoeffs = dims[2];
            jpl->nutation = Get_DoubleDataset_3D( file, "/Nutations", dims);
            jpl->nutation_nvals = dims[0]; jpl->nutation_naxes = dims[1]; jpl->nutation_ncoeffs = dims[2];
            jpl->LibrNutAlloced = TRUE;
        }
        JDparams = Get_DoubleDataset_1D( file, "/JDparams", dims);
        jpl->jalpha = JDparams[0];
        jpl->jomega = JDparams[1];
        jpl->jdelta = JDparams[2];
        LGM_ARRAY_1D_FREE( JDparams );
        //now close up
        status = H5Fclose( file );

    } else {

        printf("Problem encountered finding/opening the specified definitive ephemeris file (%s)\n", JPLephemFile);
        exit(2);

    }

}


void Lgm_JPLephem_setup_object( int objName, Lgm_JPLephemInfo *jpl, Lgm_JPLephemBundle *bundle ) {
    
    double days_per_set, offset, tdb;
    int    index, number_of_sets, number_of_axes, coefficient_count;
    int    ii, jj, omegas;
    double vx, wx, div, mod, floordiv, t1, twot1;
    double *T, ***cheby, **coefficients;

    tdb = bundle->tdb;
    if ((tdb < jpl->jalpha) || (tdb > jpl->jomega)) { 
        printf("Time (JD) is %15.8lf but must be between %15.8lf and %15.8lf\n", tdb, jpl->jalpha, jpl->jomega);
        exit(-1);
        }

    /* get DEXXX specific info */
    number_of_sets = Lgm_JPL_getNSets( objName, jpl);
    number_of_axes = Lgm_JPL_getNAxes( objName, jpl);
    coefficient_count = Lgm_JPL_getNCoeffs( objName, jpl);

    if ((number_of_sets == -1) || (number_of_axes == -1) || (coefficient_count == -1)) {
        printf("Invalid number of sets (%d), axes (%d) or coeffs (%d) for object %d\n", number_of_sets, number_of_axes, coefficient_count, objName);
        exit(-1);
        }
    days_per_set = (jpl->jomega - jpl->jalpha) / number_of_sets;
    cheby = Lgm_JPL_getCoeffSet(objName, jpl);

    /* borrow divmod implementation from Python to get index and offset */
    vx = tdb - jpl->jalpha;
    wx = days_per_set;
    mod = fmod(vx, wx);
    div = (vx - mod) / wx;
    if (mod) {
        /* ensure the remainder has the same sign as the denominator */
        if ((wx < 0) != (mod < 0)) {
            mod += wx;
            div -= 1.0;
        }
    }
    else {
        mod *= mod;  /* hide "mod = +0" from optimizer */
        if (wx < 0.0)
            mod = -mod;
    }
    /* snap quotient to nearest integral value */
    if (div) {
        floordiv = floor(div);
        if (div - floordiv > 0.5)
            floordiv += 1.0;
    }
    else {
        /* div is zero - get the same sign as the true quotient */
        div *= div;             /* hide "div = +0" from optimizers */
        floordiv = div * vx / wx; /* zero w/ sign of vx/wx */
    }
    index = (int)floordiv;
    offset = mod;

    /* extract right set of Chebyshev coefficients */
    omegas = (index==number_of_sets) ? 1 : 0;
    /* if at the end of the interval, roll back to the last set of coefficients and add offset */
    index = (omegas) ? index-1 : index;
    offset = (omegas) ? offset+days_per_set : offset;
    
    /* TODO: Check whether this can just be done with a straight assignment like
     * coefficients = cheby[index];
     */
    LGM_ARRAY_2D( coefficients, number_of_axes, coefficient_count, double );
    for (ii=0; ii<number_of_axes; ii++) {
        for (jj=0; jj<coefficient_count; jj++) {
                coefficients[ii][jj] = cheby[index][ii][jj];
            }
        }
    
    /* Chebyshev recurrence */
    LGM_ARRAY_1D( T, coefficient_count, double );
    T[0] = 1.0;
    t1 = 2.0 * offset / days_per_set - 1.0;
    twot1 = t1 + t1;
    T[1] = t1;
    for (ii=2; ii<coefficient_count; ii++) {
        T[ii] = twot1 * T[ii-1] - T[ii-2];
        }

    /* Stuff into a small struct */
    bundle->coeffs = coefficients;
    bundle->coefficient_count = coefficient_count;
    bundle->days_per_set = days_per_set;
    bundle->T = T;
    bundle->twot1 = twot1;
    bundle->TAlloced = TRUE;
    bundle->coeffsAlloced = TRUE;

    return;
    }


Lgm_JPLephemBundle *Lgm_InitJPLephemBundle( double tdb ) {
    Lgm_JPLephemBundle  *bundle = (Lgm_JPLephemBundle *) calloc (1, sizeof(Lgm_JPLephemBundle));
    bundle->tdb = tdb;
    bundle->TAlloced = FALSE;
    bundle->coeffsAlloced = FALSE;
    return bundle;
    }


void Lgm_FreeJPLephemBundle( Lgm_JPLephemBundle *bundle ) {
    if (bundle->coeffsAlloced) { LGM_ARRAY_2D_FREE( bundle->coeffs ); }
    if (bundle->TAlloced) { LGM_ARRAY_1D_FREE( bundle->T ); }
    free( bundle );
    }


void Lgm_JPLephem_position( double tdb, int objName, Lgm_JPLephemInfo *jpl, Lgm_Vector *position) {
    int i, j;
    double dum;
    Lgm_JPLephemBundle *bundle = Lgm_InitJPLephemBundle( tdb );

    //if Earth (recursive)
    if (objName==LGM_DE_EARTH) {
        Lgm_Vector EMbary, MoonGCRF;
        double earth_fac;
        Lgm_JPLephem_position( tdb, LGM_DE_EARTHMOON, jpl, &EMbary);
        Lgm_JPLephem_position( tdb, LGM_DE_MOON, jpl, &MoonGCRF);
        
        earth_fac = 1.0/(1.0+LGM_DE421_EMRAT); //TODO: fix hardcoding so that this is generic to any DEXXX loaded
        position->x = EMbary.x - MoonGCRF.x * earth_fac;
        position->y = EMbary.y - MoonGCRF.y * earth_fac;
        position->z = EMbary.z - MoonGCRF.z * earth_fac;
        }
    else if (objName==LGM_DE_MOON_ICRF) {
        Lgm_Vector EMbary, MoonGCRF;
        double moon_fac;
        Lgm_JPLephem_position( tdb, LGM_DE_EARTHMOON, jpl, &EMbary);
        Lgm_JPLephem_position( tdb, LGM_DE_MOON, jpl, &MoonGCRF);
        
        moon_fac = LGM_DE421_EMRAT/(1.0+LGM_DE421_EMRAT); //TODO: fix hardcoding so that this is generic to any DEXXX loaded
        position->x = EMbary.x + MoonGCRF.x * moon_fac;
        position->y = EMbary.y + MoonGCRF.y * moon_fac;
        position->z = EMbary.z + MoonGCRF.z * moon_fac;
        }
    else {
        //precalculate
        Lgm_JPLephem_setup_object( objName, jpl, bundle );

        //calculate positions
        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += bundle->T[i] * bundle->coeffs[0][i];
            }
        position->x = dum;

        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += bundle->T[i] * bundle->coeffs[1][i];
            }
        position->y = dum;

        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += bundle->T[i] * bundle->coeffs[2][i];
            }
        position->z = dum;

        }
    Lgm_FreeJPLephemBundle( bundle );
    }

void Lgm_JPL_getSunVector ( double tdb, Lgm_JPLephemInfo *jpl, Lgm_Vector *position ) {
    Lgm_Vector EarthICRF, SunICRF;
    // get positions of Sun & Earth
    Lgm_JPLephem_position( tdb, LGM_DE_EARTH, jpl, &EarthICRF);
    Lgm_JPLephem_position( tdb, LGM_DE_SUN, jpl, &SunICRF);

    // Calculate Earth-Sun vector 
    position->x = SunICRF.x - EarthICRF.x;
    position->y = SunICRF.y - EarthICRF.y;
    position->z = SunICRF.z - EarthICRF.z;
    }


void Lgm_JPLephem_velocity( double tdb, int objName, Lgm_JPLephemInfo *jpl, Lgm_Vector *velocity) {
    int i, j, ii;
    double *dT;
    double dum;
    Lgm_JPLephemBundle *bundle = Lgm_InitJPLephemBundle( tdb );

    //if Earth (recursive)
    if (objName==LGM_DE_EARTH) {
        Lgm_Vector EMbary, MoonGCRF;
        double earth_fac;
        Lgm_JPLephem_velocity( tdb, LGM_DE_EARTHMOON, jpl, &EMbary);
        Lgm_JPLephem_velocity( tdb, LGM_DE_MOON, jpl, &MoonGCRF);
        
        earth_fac = 1.0/(1.0+LGM_DE421_EMRAT); //TODO: fix hardcoding so that this is generic to any DEXXX loaded
        velocity->x = EMbary.x - MoonGCRF.x * earth_fac;
        velocity->y = EMbary.y - MoonGCRF.y * earth_fac;
        velocity->z = EMbary.z - MoonGCRF.z * earth_fac;
        }
    else if (objName==LGM_DE_MOON_ICRF) {
        Lgm_Vector EMbary, MoonGCRF;
        double moon_fac;
        Lgm_JPLephem_velocity( tdb, LGM_DE_EARTHMOON, jpl, &EMbary);
        Lgm_JPLephem_velocity( tdb, LGM_DE_MOON, jpl, &MoonGCRF);
        
        moon_fac = LGM_DE421_EMRAT/(1.0+LGM_DE421_EMRAT); //TODO: fix hardcoding so that this is generic to any DEXXX loaded
        velocity->x = EMbary.x + MoonGCRF.x * moon_fac;
        velocity->y = EMbary.y + MoonGCRF.y * moon_fac;
        velocity->z = EMbary.z + MoonGCRF.z * moon_fac;
        }
    else {
        //precalculate
        Lgm_JPLephem_setup_object( objName, jpl, bundle );
        /* Chebyshev derivative */
        LGM_ARRAY_1D( dT, bundle->coefficient_count, double );
        dT[0] = 0.0;
        dT[1] = 1.0;
        dT[2] = bundle->twot1 + bundle->twot1;
        for (ii=3; ii<bundle->coefficient_count; ii++) {
            dT[ii] = bundle->twot1 * dT[ii-1] - dT[ii-2] + bundle->T[ii-1] + bundle->T[ii-1];
            }
        for (ii=0; ii<bundle->coefficient_count; ii++) {
            dT[ii] *= 2.0;
            dT[ii] /= bundle->days_per_set; //puts units to km/day
            }

        //calculate velocities
        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += dT[i] * bundle->coeffs[0][i];
            }
        velocity->x = dum;

        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += dT[i] * bundle->coeffs[1][i];
            }
        velocity->y = dum;

        dum = 0;
        for (i=0; i<bundle->coefficient_count; i++) {
            dum += dT[i] * bundle->coeffs[2][i];
            }
        velocity->z = dum;

        }
    Lgm_FreeJPLephemBundle( bundle );
    }

int Lgm_JPL_getNSets( int objName, Lgm_JPLephemInfo *jpl) {
    int nvals;
    switch (objName) {
        case LGM_DE_SUN:
            nvals = (jpl->getSun) ? jpl->sun_nvals: -1;
            break;
        case LGM_DE_EARTHMOON:
            nvals = (jpl->getEarth) ? jpl->earthmoon_nvals: -1;
            break;
        case LGM_DE_MOON:
            nvals = (jpl->getEarth) ? jpl->moon_wrt_earth_nvals: -1;
            break;
        case LGM_DE_MERCURY:
            nvals = (jpl->getInnerPlanets) ? jpl->mercury_nvals: -1;
            break;
        case LGM_DE_VENUS:
            nvals = (jpl->getInnerPlanets) ? jpl->venus_nvals: -1;
            break;
        case LGM_DE_MARS:
            nvals = (jpl->getInnerPlanets) ? jpl->mars_nvals: -1;
            break;
        case LGM_DE_JUPITER:
            nvals = (jpl->getOuterPlanets) ? jpl->jupiter_nvals: -1;
            break;
        case LGM_DE_SATURN:
            nvals = (jpl->getOuterPlanets) ? jpl->saturn_nvals: -1;
            break;
        case LGM_DE_URANUS:
            nvals = (jpl->getOuterPlanets) ? jpl->uranus_nvals: -1;
            break;
        case LGM_DE_NEPTUNE:
            nvals = (jpl->getOuterPlanets) ? jpl->neptune_nvals: -1;
            break;
        case LGM_DE_PLUTO:
            nvals = (jpl->getOuterPlanets) ? jpl->pluto_nvals: -1;
            break;
        case LGM_DE_LIBRATION:
            nvals = (jpl->getLibrationNutation) ? jpl->libration_nvals: -1;
            break;
        case LGM_DE_NUTATION:
            nvals = (jpl->getLibrationNutation) ? jpl->nutation_nvals: -1;
            break;
        }
    return nvals;
    }

int Lgm_JPL_getNAxes( int objName, Lgm_JPLephemInfo *jpl) {
    int naxes;
    switch (objName) {
        case LGM_DE_SUN:
            naxes = (jpl->getSun) ? jpl->sun_naxes: -1;
            break;
        case LGM_DE_EARTHMOON:
            naxes = (jpl->getEarth) ? jpl->earthmoon_naxes: -1;
            break;
        case LGM_DE_MOON:
            naxes = (jpl->getEarth) ? jpl->moon_wrt_earth_naxes: -1;
            break;
        case LGM_DE_MERCURY:
            naxes = (jpl->getInnerPlanets) ? jpl->mercury_naxes: -1;
            break;
        case LGM_DE_VENUS:
            naxes = (jpl->getInnerPlanets) ? jpl->venus_naxes: -1;
            break;
        case LGM_DE_MARS:
            naxes = (jpl->getInnerPlanets) ? jpl->mars_naxes: -1;
            break;
        case LGM_DE_JUPITER:
            naxes = (jpl->getOuterPlanets) ? jpl->jupiter_naxes: -1;
            break;
        case LGM_DE_SATURN:
            naxes = (jpl->getOuterPlanets) ? jpl->saturn_naxes: -1;
            break;
        case LGM_DE_URANUS:
            naxes = (jpl->getOuterPlanets) ? jpl->uranus_naxes: -1;
            break;
        case LGM_DE_NEPTUNE:
            naxes = (jpl->getOuterPlanets) ? jpl->neptune_naxes: -1;
            break;
        case LGM_DE_PLUTO:
            naxes = (jpl->getOuterPlanets) ? jpl->pluto_naxes: -1;
            break;
        case LGM_DE_LIBRATION:
            naxes = (jpl->getLibrationNutation) ? jpl->libration_naxes: -1;
            break;
        case LGM_DE_NUTATION:
            naxes = (jpl->getLibrationNutation) ? jpl->nutation_naxes: -1;
            break;
        default:
            naxes = -1 ;
            break;
        }
    return naxes;
    }


int Lgm_JPL_getNCoeffs( int objName, Lgm_JPLephemInfo *jpl) {
    int ncoeffs;
    switch (objName) {
        case LGM_DE_SUN:
            ncoeffs = (jpl->getSun) ? jpl->sun_ncoeffs: -1;
            break;
        case LGM_DE_EARTHMOON:
            ncoeffs = (jpl->getEarth) ? jpl->earthmoon_ncoeffs: -1;
            break;
        case LGM_DE_MOON:
            ncoeffs = (jpl->getEarth) ? jpl->moon_wrt_earth_ncoeffs: -1;
            break;
        case LGM_DE_MERCURY:
            ncoeffs = (jpl->getInnerPlanets) ? jpl->mercury_ncoeffs: -1;
            break;
        case LGM_DE_VENUS:
            ncoeffs = (jpl->getInnerPlanets) ? jpl->venus_ncoeffs: -1;
            break;
        case LGM_DE_MARS:
            ncoeffs = (jpl->getInnerPlanets) ? jpl->mars_ncoeffs: -1;
            break;
        case LGM_DE_JUPITER:
            ncoeffs = (jpl->getOuterPlanets) ? jpl->jupiter_ncoeffs: -1;
            break;
        case LGM_DE_SATURN:
            ncoeffs = (jpl->getOuterPlanets) ? jpl->saturn_ncoeffs: -1;
            break;
        case LGM_DE_URANUS:
            ncoeffs = (jpl->getOuterPlanets) ? jpl->uranus_ncoeffs: -1;
            break;
        case LGM_DE_NEPTUNE:
            ncoeffs = (jpl->getOuterPlanets) ? jpl->neptune_ncoeffs: -1;
            break;
        case LGM_DE_PLUTO:
            ncoeffs = (jpl->getOuterPlanets) ? jpl->pluto_ncoeffs: -1;
            break;
        case LGM_DE_LIBRATION:
            ncoeffs = (jpl->getLibrationNutation) ? jpl->libration_ncoeffs: -1;
            break;
        case LGM_DE_NUTATION:
            ncoeffs = (jpl->getLibrationNutation) ? jpl->nutation_ncoeffs: -1;
            break;
        default:
            ncoeffs = -1 ;
            break;
        }
    return ncoeffs;
    }

double ***Lgm_JPL_getCoeffSet(int objName, Lgm_JPLephemInfo *jpl){
    double ***coeffset;
    switch (objName) {
        case LGM_DE_SUN:
            coeffset = (jpl->getSun) ? jpl->sun: NULL;
            break;
        case LGM_DE_EARTHMOON:
            coeffset = (jpl->getEarth) ? jpl->earthmoon: NULL;
            break;
        case LGM_DE_MOON:
            coeffset = (jpl->getEarth) ? jpl->moon_wrt_earth: NULL;
            break;
        case LGM_DE_MERCURY:
            coeffset = (jpl->getInnerPlanets) ? jpl->mercury: NULL;
            break;
        case LGM_DE_VENUS:
            coeffset = (jpl->getInnerPlanets) ? jpl->venus: NULL;
            break;
        case LGM_DE_MARS:
            coeffset = (jpl->getInnerPlanets) ? jpl->mars: NULL;
            break;
        case LGM_DE_JUPITER:
            coeffset = (jpl->getOuterPlanets) ? jpl->jupiter: NULL;
            break;
        case LGM_DE_SATURN:
            coeffset = (jpl->getOuterPlanets) ? jpl->saturn: NULL;
            break;
        case LGM_DE_URANUS:
            coeffset = (jpl->getOuterPlanets) ? jpl->uranus: NULL;
            break;
        case LGM_DE_NEPTUNE:
            coeffset = (jpl->getOuterPlanets) ? jpl->neptune: NULL;
            break;
        case LGM_DE_PLUTO:
            coeffset = (jpl->getOuterPlanets) ? jpl->pluto: NULL;
            break;
        case LGM_DE_LIBRATION:
            coeffset = (jpl->getLibrationNutation) ? jpl->libration: NULL;
            break;
        case LGM_DE_NUTATION:
            coeffset = (jpl->getLibrationNutation) ? jpl->nutation: NULL;
            break;
        default:
            coeffset = NULL ;
            break;
        }
    return coeffset;
    }
