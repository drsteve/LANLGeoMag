#include "Lgm/Lgm_MagModelInfo.h"


/*
 *  A convenience routine to set some default values
 */
/*
 * BAL 19Jan2011 split into 2 sections for python wrapping without the calloc
 */
Lgm_MagModelInfo *Lgm_InitMagInfo( ) {
    Lgm_MagModelInfo  *MagInfo;
    MagInfo = (Lgm_MagModelInfo *) calloc (1, sizeof(*MagInfo));
    Lgm_InitMagInfoDefaults(MagInfo);
    return(MagInfo);
}


void Lgm_InitMagInfoDefaults( Lgm_MagModelInfo  *MagInfo ) {

    MagInfo->Bfield = Lgm_B_T89;
    MagInfo->InternalModel = LGM_IGRF;

    MagInfo->c     = Lgm_init_ctrans( 0 );
    MagInfo->Kp    = 5;
    MagInfo->P     = 2.1; // SW pressure in nPa
    MagInfo->nFunc = 0;

    MagInfo->B0    = 1.00;   // Should nominally be 1.0 See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B1    = 0.8100; // See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B2    = 0.4065; // See page 30 of Schultz and Lanzerotti, [1974]

    MagInfo->SavePoints = 0;
    MagInfo->Hmax       = 1.0;

    MagInfo->UseInterpRoutines = TRUE;
    Lgm_Set_Open_Limits( MagInfo, -80.0, 30.0, -40.0, 40.0, -40.0, 40.0 );

    MagInfo->Lgm_I_integrand_JumpMethod = LGM_ABSOLUTE_JUMP_METHOD;

    /*
     *  Some inits for MagStep
     */
    MagInfo->Lgm_MagStep_FirstTimeThrough = TRUE;
    MagInfo->Lgm_MagStep_eps_old          = -1.0;

//    gsl_set_error_handler_off(); // Turn off gsl default error handler

    /*
     * Set some default tolerances
     */
    MagInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
    MagInfo->Lgm_MagFlux_Integrator_epsrel = 1e-5;

    MagInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
    MagInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

    MagInfo->Lgm_I_Integrator_epsrel = 0.0;
    MagInfo->Lgm_I_Integrator_epsabs = 1e-3;
    MagInfo->Lgm_I_Integrator = DQAGS;

    MagInfo->Lgm_Sb_Integrator_epsrel = 0.0;
    MagInfo->Lgm_Sb_Integrator_epsabs = 1e-3;
    MagInfo->Lgm_Sb_Integrator = DQAGP; // not changeable (yet...)

    MagInfo->Lgm_FindBmRadius_Tol = 1e-10;
    MagInfo->Lgm_FindShellLine_I_Tol = 1e-3;
    MagInfo->Lgm_TraceToMirrorPoint_Tol = 1e-7;

    /*
     * Bounce Loss Cone Height
     */
    MagInfo->Lgm_LossConeHeight = 100.0; // km above the Earth Ellipsoid.

    /*
     * Inits for OP77
     */
    MagInfo->OP77_TILTL = 99.0;
}



void Lgm_FreeMagInfo( Lgm_MagModelInfo  *Info ) {

    Lgm_free_ctrans( Info->c );
    free( Info );

}




/*
 *  The Lgm_MagModelInfo structure has pointers in it, so simple
 *  asignments (e.g. *t = *s) are dangerous. Here we make sure that
 *  the target gets an independent copy of the structure.
 */
Lgm_MagModelInfo *Lgm_CopyMagInfo( Lgm_MagModelInfo *s ) {

    Lgm_MagModelInfo *t;


    if ( s == NULL) {
        printf("Lgm_CopyMagInfo: Error, source structure is NULL\n");
        return;
    }

    t = (Lgm_MagModelInfo *)calloc( 1, sizeof(Lgm_MagModelInfo) );

    // memcpy's args are (dest, src, size)
    memcpy( t, s, sizeof(Lgm_MagModelInfo) );


    /*
     *  Now, copy the Lgm_CTrans struct properly.
     */
    t->c = Lgm_CopyCTrans( s->c );


    /*
     *  Lets also assume that GSL interp stuff should be NULL
     */
    //t->acc    = (gsl_interp_accel *)NULL;
    //t->spline = (gsl_spline *)NULL;

    return( t );

}



/*
 * Testing...
 * This type of model for setting things is more in the style of OO-programming
 * Could return error if pointer is NULL???
 */
void Lgm_MagModelInfo_Set_Psw( double Psw, Lgm_MagModelInfo *m ) {
    m->P = Psw;
}
void Lgm_MagModelInfo_Set_Kp( double Kp, Lgm_MagModelInfo *m ) {
    m->Kp = Kp;
}
void Lgm_MagModelInfo_Set_MagModel( Lgm_MagModelInfo *m, int InternalModel, int ExternalModel ){

    m->InternalModel = InternalModel;
    m->ExternalModel = ExternalModel;
    
    switch ( m->ExternalModel ) {

        case LGM_EXTMODEL_NULL:
                                /*
                                 * This is a bit kludgey. If NULL is given
                                 * assume the user wants only an internal
                                 * field. I.e. set m->Bfield to be equal to the
                                 * Internal model given.
                                 */
                                if ( InternalModel == LGM_CDIP ) {
                                    m->Bfield = Lgm_B_cdip;
                                } else if ( InternalModel == LGM_EDIP ) {
                                    m->Bfield = Lgm_B_edip;
                                } else {
                                    m->Bfield = Lgm_B_igrf;
                                }
                                break;

        case LGM_EXTMODEL_T87:
                                m->Bfield = Lgm_B_T87;
                                break;

        case LGM_EXTMODEL_T89:
                                m->Bfield = Lgm_B_T89;
                                break;

/*
        case LGM_EXTMODEL_T96:
                                m->Bfield = Lgm_B_T96;
                                break;
*/

        case LGM_EXTMODEL_T01S:
                                m->Bfield = Lgm_B_T01S;
                                break;

        case LGM_EXTMODEL_TS04:
                                m->Bfield = Lgm_B_TS04;
                                break;

        case LGM_EXTMODEL_OP77:
                                m->Bfield = Lgm_B_OP77;
                                break;

    }

}

/*
 * Some setters
 */
void Lgm_Set_Octree_kNN_InterpMethod( Lgm_MagModelInfo *m, int Method ) {
    m->Octree_kNN_InterpMethod = Method;
    return;
}
void Lgm_Set_Octree_kNN_k( Lgm_MagModelInfo *m, int k ) {
    m->Octree_kNN_k = k;
    return;
}
void Lgm_Set_Octree_kNN_MaxDist( Lgm_MagModelInfo *m, double MaxDist ) {
    m->Octree_kNN_MaxDist = MaxDist; // This maxdist is in physical units
    return;
}


void Lgm_Set_Open_Limits( Lgm_MagModelInfo *m, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax ){
    m->OpenLimit_xMin = xmin;
    m->OpenLimit_xMax = xmax;
    m->OpenLimit_yMin = ymin;
    m->OpenLimit_yMax = ymax;
    m->OpenLimit_zMin = zmin;
    m->OpenLimit_zMax = zmax;
    return;
}


void Lgm_Set_LossConeHeight( Lgm_MagModelInfo *m, double LossConeHeight ) {
    m->Lgm_LossConeHeight = LossConeHeight;
    return;
}

void Lgm_Set_Lgm_B_igrf(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_igrf;
}

void Lgm_Set_Lgm_B_T01S(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_T01S;
}

void Lgm_Set_gm_B_TS04(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_TS04;
}

void Lgm_Set_Lgm_B_T89(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_T89;
}

void Lgm_Set_Lgm_B_OP77(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_OP77;
}





/*
 *   $Id: Lgm_InitMagInfo.c 139 2011-01-27 21:01:34Z mgh $
 */
