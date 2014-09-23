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

    MagInfo->AllocedSplines = FALSE;

    MagInfo->Bfield = Lgm_B_T89;
    MagInfo->InternalModel = LGM_IGRF;

    MagInfo->c     = Lgm_init_ctrans( 0 );
    MagInfo->fKp   = 2.0;
    MagInfo->Kp    = 2;
    MagInfo->P     = 2.1; // SW pressure in nPa
    MagInfo->nFunc = 0;
    MagInfo->Dst   = -5.0;

    MagInfo->B0    = 1.00;   // Should nominally be 1.0 See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B1    = 0.8100; // See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B2    = 0.4065; // See page 30 of Schultz and Lanzerotti, [1974]

    MagInfo->SavePoints = 0;
    MagInfo->Hmax       = 1.0;

    MagInfo->ComputeSb0 = FALSE;

    MagInfo->UseInterpRoutines = TRUE;
    Lgm_Set_Open_Limits( MagInfo, -80.0, 30.0, -40.0, 40.0, -40.0, 40.0 );

    MagInfo->Lgm_I_integrand_JumpMethod = LGM_ABSOLUTE_JUMP_METHOD;

    /*
     * Default to Bulirsch-Stoer ODE method for FL tracing (i.e. LGM_MAGSTEP_ODE_BS)
     * Could also use LGM_MAGSTEP_ODE_RK5
     */
    MagInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;

    /*
     *  Some inits for MagStep_BS
     */
    MagInfo->Lgm_MagStep_BS_FirstTimeThrough = TRUE;
    MagInfo->Lgm_MagStep_BS_eps_old          = -1.0;
    MagInfo->Lgm_MagStep_BS_Eps              = 1e-7;

    /*
     *  Some inits for new MagStep_BS
     */
    MagInfo->Lgm_MagStep_BS_first_step  = TRUE;
    MagInfo->Lgm_MagStep_BS_last_step   = FALSE;
    MagInfo->Lgm_MagStep_BS_reject      = FALSE;
    MagInfo->Lgm_MagStep_BS_prev_reject = FALSE;
    MagInfo->Lgm_MagStep_BS_atol        = 1e-5;
    MagInfo->Lgm_MagStep_BS_rtol        = 1e-5;

    /*
     *  Some inits for MagStep_RK5
     */
    MagInfo->Lgm_MagStep_RK5_FirstTimeThrough = TRUE;
    MagInfo->Lgm_MagStep_RK5_snew             =  0.0;
    MagInfo->Lgm_MagStep_RK5_MaxCount         =  50;
    MagInfo->Lgm_MagStep_RK5_Safety           =  0.9;
    MagInfo->Lgm_MagStep_RK5_pGrow            = -0.2;
    MagInfo->Lgm_MagStep_RK5_pShrnk           = -0.25;
    MagInfo->Lgm_MagStep_RK5_ErrCon           = pow( 5.0/MagInfo->Lgm_MagStep_RK5_Safety, 1.0/MagInfo->Lgm_MagStep_RK5_pGrow);
    MagInfo->Lgm_MagStep_RK5_Eps              = 1e-5;

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
    MagInfo->Lgm_Sb_Integrator_epsabs = 1e-4;
    MagInfo->Lgm_Sb_Integrator = DQAGP; // not changeable (yet...)

    MagInfo->Lgm_FindBmRadius_Tol = 1e-10;
    MagInfo->Lgm_FindShellLine_I_Tol = 1e-3;
    MagInfo->Lgm_TraceToMirrorPoint_Tol = 1e-7;
    MagInfo->Lgm_TraceToEarth_Tol = 1e-7;
    MagInfo->Lgm_TraceToBmin_Tol  = 1e-7;
    MagInfo->Lgm_TraceLine_Tol    = 1e-7;



    /*
     * Step sizes for pre-tracing routines (e.g. TraceLine() type routines that
     * save quantities into arrays for later use in interp schemes).
     */
    MagInfo->MaxDiv = 0.1;
    MagInfo->nDivs  = 200;

    /*
     * Bounce Loss Cone Height
     */
    MagInfo->Lgm_LossConeHeight = 100.0; // km above the Earth Ellipsoid.

    /*
     * Inits for OP77
     */
    MagInfo->OP77_TILTL = 99.0;


    /*
     * Inits for T96
     */
    Lgm_Init_T96( &MagInfo->T96_Info );

    /*
     * Inits for T01
     */
    Lgm_Init_T01S( &MagInfo->T01_Info );

    /*
     * Inits for TS04
     */
    Lgm_Init_TS04( &MagInfo->TS04_Info );


    /*
     * Inits for TS07
     */
    Lgm_Init_TS07( &MagInfo->TS07_Info );


    /*
     *  Inits for Octree stuff
     */
    //MagInfo->Octree_kNN_InterpMethod = XXX; // need to redo this
    MagInfo->Octree_kNN_k = 4;
    MagInfo->Octree_kNN_MaxDist2 = 1.0*1.0;
    MagInfo->Octree_kNN_Alloced = 0;


    /*
     *  Initialize hash table used in Lgm_B_FromScatteredData*()
     */
    MagInfo->rbf_ht         = NULL;
    MagInfo->rbf_ht_alloced = FALSE;
    MagInfo->RBF_nHashFinds = 0;
    MagInfo->RBF_nHashAdds  = 0;


    /*
     * Allocate space to FastPow structure
     * And initialize it...
     */
//    MagInfo->f = Lgm_InitFastPow();
//    powFastSetTable( MagInfo->f );


}


void Lgm_FreeMagInfo_children( Lgm_MagModelInfo  *Info ) {

    Lgm_DeAllocate_TS07( &(Info->TS07_Info) );
    Lgm_free_ctrans( Info->c );



//    Lgm_FreeFastPow( Info->f );

    //what about any splines that may have been alloc'd?


}


void Lgm_FreeMagInfo( Lgm_MagModelInfo  *Info ) {
    Lgm_FreeMagInfo_children( Info );
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
        return((Lgm_MagModelInfo *) NULL);
    }

    t = (Lgm_MagModelInfo *)calloc( 1, sizeof(Lgm_MagModelInfo) );

    // memcpy's args are (dest, src, size)
    memcpy( t, s, sizeof(Lgm_MagModelInfo) );


    /*
     *  t->c and s->c now point at the same thing. We dont want that.
     *  Now, copy the Lgm_CTrans struct properly.
     */
    t->c = Lgm_CopyCTrans( s->c );




    /*
     *  Lets also assume that GSL interp stuff should be NULL
     *  This routine wont yet properly copy over the spline stuff.
     *  Need to be careful if you want an independent copy of that stuff...
     *  BE careful free these things in copies -- it'll free them for all other copies too.
     */
    //t->acc    = (gsl_interp_accel *)NULL;
    //t->spline = (gsl_spline *)NULL;

    // octree stuff is also not copied correctly...

    if ( s->Octree_Alloced ) {
        //t->Octree         = Lgm_CopyOctree( s->Octree );
        t->Octree         = s->Octree;
        t->Octree_Alloced = TRUE;
    } else {
        t->Octree         = NULL;
        t->Octree_Alloced = FALSE;
    }



    // TEMP KLUDGE
    t->TS07_Info.ArraysAlloced = FALSE;


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
void Lgm_MagModelInfo_Set_MagModel( int InternalModel, int ExternalModel, Lgm_MagModelInfo *m ){

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
                                } else if ( InternalModel == LGM_DUNGEY ) {
                                    m->Bfield = Lgm_B_Dungey;
                                } else {
                                    m->Bfield = Lgm_B_igrf;
                                }
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_TU82:
                                m->Bfield = Lgm_B_TU82;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;
        case LGM_EXTMODEL_OP88:
                                m->Bfield = Lgm_B_OP88;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;
        case LGM_EXTMODEL_T87:
                                m->Bfield = Lgm_B_T87;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_T89:
                                m->Bfield = Lgm_B_T89;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_T89c:
                                m->Bfield = Lgm_B_T89c;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_T96:
                                m->Bfield = Lgm_B_T96;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_T02:
                                m->Bfield = Lgm_B_T02;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_T01S:
                                m->Bfield = Lgm_B_T01S;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_TS04:
                                m->Bfield = Lgm_B_TS04;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;
        case LGM_EXTMODEL_TS07:
                                m->Bfield = Lgm_B_TS07;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_OP77:
                                m->Bfield = Lgm_B_OP77;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_SCATTERED_DATA:
                                m->Bfield = Lgm_B_FromScatteredData;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        case LGM_EXTMODEL_SCATTERED_DATA2:
                                m->Bfield = Lgm_B_FromScatteredData2;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                break;

        default:
                                printf("No such B field model\n");
                                exit(-1);
                                break;


    }

}

/*
 * Some setters
 */
void Lgm_MagModelInfo_Set_Octree( Lgm_Octree *Octree, int k, Lgm_MagModelInfo *m ) {
    m->Octree         = Octree;
    m->Octree_Alloced = TRUE;
    Lgm_Set_Octree_kNN_k( m, k );
    return;
}
void Lgm_Set_Octree_kNN_InterpMethod( Lgm_MagModelInfo *m, int Method ) {
    m->Octree_kNN_InterpMethod = Method;
    return;
}
void Lgm_Set_Octree_kNN_k( Lgm_MagModelInfo *m, int k ) {
    m->Octree_kNN_k = k;
    return;
}
void Lgm_Set_Octree_kNN_MaxDist2( Lgm_MagModelInfo *m, double MaxDist2 ) {
    m->Octree_kNN_MaxDist2 = MaxDist2; // This maxdist is in physical units
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

void Lgm_Set_Lgm_B_TS04(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_TS04;
}

void Lgm_Set_Lgm_B_TS07(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_TS07;
}

void Lgm_Set_Lgm_B_T96(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_T96;
}

void Lgm_Set_Lgm_B_T89(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_T89;
}

void Lgm_Set_Lgm_B_T89c(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_T89c;
}

void Lgm_Set_Lgm_B_OP77(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_OP77;
}

void Lgm_Set_Lgm_B_Dungey(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_DUNGEY;
    MagInfo->Bfield = Lgm_B_Dungey;
}

void Lgm_Set_Lgm_B_cdip(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_CDIP;
    MagInfo->Bfield = Lgm_B_cdip;
}

void Lgm_Set_Lgm_B_edip(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_EDIP;
    MagInfo->Bfield = Lgm_B_edip;
}

void Lgm_Set_Lgm_B_cdip_InternalModel(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_CDIP;
}

void Lgm_Set_Lgm_B_edip_InternalModel(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_EDIP;
}

void Lgm_Set_Lgm_B_IGRF_InternalModel(Lgm_MagModelInfo *MagInfo) {
    MagInfo->InternalModel = LGM_IGRF;
}

