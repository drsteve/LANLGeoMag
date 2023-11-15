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

    MagInfo->Bfield = Lgm_B_T89c;
    MagInfo->InternalModel = LGM_IGRF;

    MagInfo->c     = Lgm_init_ctrans( 0 );
    MagInfo->nFunc = 0;

    MagInfo->B0    = 1.00;   // Should nominally be 1.0 See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B1    = 0.8100; // See page 30 of Schultz and Lanzerotti, [1974]
    MagInfo->B2    = 0.4065; // See page 30 of Schultz and Lanzerotti, [1974]



    /*
     * QinDenton Defaults
     * Make sure we have at least something for default values here.
     */
    MagInfo->Dst   =  -5.0; // nT
    MagInfo->fKp   =   2.0; // dimensionless
    MagInfo->Kp    =     2; // dimensionless
    MagInfo->P     =   2.1; // Pressure in nPa
    MagInfo->By    =   5.0; // IMF By in nT why are there two of these? (Maybe QinDenton codes introduced the new set?)
    MagInfo->Bz    =  -5.0; // IMF Bz in nT
    MagInfo->W[0]  =  0.44; // W1 these are avg values at status flag 2 (Table 3, Qin et al.)
    MagInfo->W[1]  =  0.42; // W2
    MagInfo->W[2]  =  0.66; // W3
    MagInfo->W[3]  =  0.48; // W4
    MagInfo->W[4]  =  0.49; // W5
    MagInfo->W[5]  =  0.91; // W6
    MagInfo->G1    =   6.0; // G1
    MagInfo->G2    =  10.0; // G2
    MagInfo->G3    =  60.0; // G3



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
    MagInfo->Lgm_MagStep_BS_rtol        = 0.0;

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

    MagInfo->Lgm_n_V_integrand_Calls = 0;
    MagInfo->Lgm_V_Integrator_epsrel = 0.0;
    MagInfo->Lgm_V_Integrator_epsabs = 1e-3;
    MagInfo->Lgm_V_Integrator = DQAGS;



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
    // Do this in the first call that sets coeff files.
    // I.e. do it in first call to Lgm_SetCoeffs_TS07()
    // Lgm_Init_TS07( &MagInfo->TS07_Info );


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
    MagInfo->rbf_ht              = NULL;
    MagInfo->rbf_ht_alloced      = FALSE;
    MagInfo->RBF_nHashFinds      = 0;
    MagInfo->RBF_nHashAdds       = 0;
    MagInfo->RBF_CompGradAndCurl = FALSE;
    MagInfo->RBF_Type            = LGM_RBF_MULTIQUADRIC;
    MagInfo->RBF_Eps             = 1.0/(4.0*4.0);
    
    MagInfo->KdTree              = NULL;
    MagInfo->KdTree_Alloced      = FALSE;
    MagInfo->KdTree_kNN_InterpMethod = 0;
    MagInfo->KdTree_kNN_k        = 12;
    MagInfo->KdTree_kNN_MaxDist2 = 1e6;

    MagInfo->KdTree_kNN          = NULL;
    MagInfo->KdTree_kNN_Alloced  = FALSE;

    MagInfo->KdTreeCopy          = FALSE;


    /*
     * Allocate space to FastPow structure
     * And initialize it...
     */
//    MagInfo->f = Lgm_InitFastPow();
//    powFastSetTable( MagInfo->f );


}


void Lgm_FreeMagInfo_children( Lgm_MagModelInfo  *Info ) {

    Lgm_DeAllocate_TS07( &(Info->TS07_Info) );
    Lgm_DeAllocate_TA16( &(Info->TA16_Info) );
    Lgm_free_ctrans( Info->c );

    if ( Info->KdTree_Alloced ) {

        if ( Info->KdTreeCopy ) {
            Lgm_KdTree_FreeLite( Info->KdTree );
        } else {
            Lgm_KdTree_Free( Info->KdTree );
        }
        Info->KdTree_Alloced = FALSE;

    }

    if ( Info->KdTree_kNN_Alloced ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
    }


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



    /* Do a "CopyLite" on the KdTree if its alloced -- this make a new independent copy with the exception of the tree data. The tree data is 
     * not copied, we just copy the pointer. To guard against freein this, we flag that this is a copy.
     */
//printf("s->KdTree_Alloced = %d\n", s->KdTree_Alloced);
//printf("s->KdTree = %s\n", s->KdTree);
    if ( s->KdTree_Alloced ) {
        t->KdTree         = Lgm_KdTree_CopyLite( s->KdTree );
        t->KdTree_Alloced = TRUE;
        t->KdTreeCopy     = TRUE;
    } else {
        t->KdTree         = NULL;
        t->KdTree_Alloced = FALSE;
        t->KdTreeCopy     = FALSE;
    }

    // dont copy the pointerj here.
    t->KdTree_kNN         = NULL;
    t->KdTree_kNN_Alloced = 0;

    // Make sure the RBF hash tables are reset.
    t->rbf_ht_alloced = FALSE;



    /*
     * Copy the TS07_Info structure
     */
    if ( s->TS07_Info.ArraysAlloced ) {
        Lgm_Copy_TS07_Info( &(t->TS07_Info), &(s->TS07_Info) );
    }



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



/*
 *  Just Calls Lgm_MagModelInfo_Set_MagModel()
 */
void Lgm_Set_MagModel( int InternalModel, int ExternalModel, Lgm_MagModelInfo *m ){
    Lgm_MagModelInfo_Set_MagModel( InternalModel, ExternalModel, m );
}

void Lgm_MagModelInfo_Set_MagModel( int InternalModel, int ExternalModel, Lgm_MagModelInfo *m ){

    m->InternalModel = InternalModel;
    m->ExternalModel = ExternalModel;

    switch( m->InternalModel ) {

        case LGM_CDIP:
                        strcpy( m->IntMagModelStr1, "CDIP" );
                        strcpy( m->IntMagModelStr2, "Centered Dipole Model." );
                        strcpy( m->IntMagModelStr3, "Reference: International Geomagnetic Reference Field: the 12th generation, Thébault et al. Earth, Planets and Space 2015, 67:79 (27 May 2015)." );
                        strcpy( m->IntMagModelStr4, "Comments: Use first 3 time-dependent IGRF coefficients." );
                        break;

        case LGM_EDIP:
                        strcpy( m->IntMagModelStr1, "EDIP" );
                        strcpy( m->IntMagModelStr2, "Eccentric Dipole Model." );
                        strcpy( m->IntMagModelStr3, "References: (1) International Geomagnetic Reference Field: the 12th generation, Thébault et al. Earth, Planets and Space 2015, 67:79 (27 May 2015); (2) Eccentric dipole coordinates for Magsat data presentation and analysis of ecternal current effects (1982), Geophys. Res. Lett., 9(4), pp. 353-356; (3) Fraser-Smith, A. C., Centered and Eccentric Geomagnetic Dipoles and Their Poles, 1600-1985, Rev. Geophys., 25, 1, pp. 1-16, 1987." );
                        strcpy( m->IntMagModelStr4, "Comments: Use first 8 time-dependent IGRF coefficients. " );
                        break;

        case LGM_DUNGEY:
                        strcpy( m->IntMagModelStr1, "DUNGEY" );
                        strcpy( m->IntMagModelStr2, "Dungey Field Model" );
                        strcpy( m->IntMagModelStr3, "Reference: See Schulz, M. (1997), Direct influence of ring current on auroral oval diameter, J. Geophys. Res., 102(A7), 14149–14154, doi:10.1029/97JA00827." );
                        strcpy( m->IntMagModelStr4, "Comments: Centered dipole + a uniform -Bz component. Not very realistic, but L* can be computed analytically, so it makes a good test for L* numerics." );
                        break;

        case LGM_JENSENCAIN1960:
                        strcpy( m->IntMagModelStr1, "DUNGEY" );
                        strcpy( m->IntMagModelStr2, "Dungey Field Model" );
                        strcpy( m->IntMagModelStr3, "Reference: See Jensen DC, and Cain, JC (1962), An interim geomagnetic field, J. Geophys. Res., 67, 3568–3569" );
                        strcpy( m->IntMagModelStr4, "Comments: Low-order approximation to IGRF valid for 1960. Used to bin some trapped particle measurements after Starfish." );
                        break;

        case LGM_IGRF:
        default:
                        strcpy( m->IntMagModelStr1, "IGRF12" );
                        strcpy( m->IntMagModelStr2, "International Geomagnetic Reference Field: the 12th generation" );
                        strcpy( m->IntMagModelStr3, "Reference: International Geomagnetic Reference Field: the 12th generation, Thébault et al. Earth, Planets and Space 2015, 67:79 (27 May 2015)." );
                        strcpy( m->IntMagModelStr4, "Comments: Uses all available coefficients for all radial distances." );
                        break;

    }
        




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
                                } else if ( InternalModel == LGM_JENSENCAIN1960 ) {
                                    m->Bfield = Lgm_B_JensenCain1960;	
                                } else {
                                    m->Bfield = Lgm_B_igrf;
                                }
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "NULL" );
                                strcpy( m->ExtMagModelStr2, "No External model - internal only" );
                                strcpy( m->ExtMagModelStr3, "Reference: N/A" );
                                strcpy( m->ExtMagModelStr4, "Comments: Running with only internal model." );
                                break;

        case LGM_EXTMODEL_TU82:
                                m->Bfield = Lgm_B_TU82;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "TU82" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko Usmanov 1982 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A. and A. V. Usmanov, Determination of the magnetospheric current system parameters and development of experimental geomagnetic field models based on data from IMP and HEOS satellites, Planet. Space Sci., 30 (1982), pp.  985-998." );
                                strcpy( m->ExtMagModelStr4, "Comments: " );

                                break;
        case LGM_EXTMODEL_OP88:
                                m->Bfield = Lgm_B_OP88;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "OP88" );
                                strcpy( m->ExtMagModelStr2, "Olson Pfitzer 1988 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Pfitzer, K. A., W. P. Olson, and T. Mogstad (1988), A time-dependent, source-driven magnetospheric magnetic field model, Eos Trans. AGU, 69(16), 426." );
                                strcpy( m->ExtMagModelStr4, "Comments: " );

                                break;
        case LGM_EXTMODEL_T87:
                                m->Bfield = Lgm_B_T87;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T87" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 1987 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A., Global quantitative models of the geomagnetic field in the cislunar magnetosphere for different disturbance levels, Planet. Space Sci., 35 (1987), pp. 1347-1358. doi:10.1016/0032-0633(87)90046-8." );
                                strcpy( m->ExtMagModelStr4, "Comments: " );
                                break;

//should we chanmge this to T89orig or some such thing to avoid confusion????
        case LGM_EXTMODEL_T89:
                                m->Bfield = Lgm_B_T89;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T89" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 1989 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A. (1989), A magnetospheric magnetic field model with a warped tail current sheet, Planet. Space Sci., 37, 5-20, doi:10.1016/0032-0633(89)90066-4.");
                                strcpy( m->ExtMagModelStr4, "Comments: The original published version of the code (with corrections to coefficients.)" );
                                break;

        case LGM_EXTMODEL_T89c:
                                m->Bfield = Lgm_B_T89c;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T89c" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 1989c Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A. (1989), A magnetospheric magnetic field model with a warped tail current sheet, Planet. Space Sci., 37, 5-20, doi:10.1016/0032-0633(89)90066-4.");
                                strcpy( m->ExtMagModelStr4, "Comments: A modified version of the originally published model. This is the version that is now most commonly referred to as T89.");
                                break;

        case LGM_EXTMODEL_T96:
                                m->Bfield = Lgm_B_T96;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T96" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 1996 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A. (1995), Modeling the Earth's magnetospheric magnetic field confined within a realistic magnetopause, J. Geophys. Res., 100(A4), 5599–5612, doi:10.1029/94JA03193." );
                                strcpy( m->ExtMagModelStr4, "Comments: ");
                                break;

        case LGM_EXTMODEL_T02:
                                m->Bfield = Lgm_B_T02;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T02" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 2002 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "References: (1) Tsyganenko, N. A., A model of the magnetosphere with a dawn-dusk asymmetry, 1, Mathematical structure, J. Geophys. Res., 107(A8), doi:10.1029/2001JA000219, 2002; (2) Tsyganenko, N. A., A model of the near magnetosphere with a dawn-dusk asymmetry, 2, Parameterization and fitting to observations, J. Geophys. Res., 107(A8), doi:10.1029/2001JA000220, 2002.");
                                strcpy( m->ExtMagModelStr4, "Comments: Also known as T01_01.");
                                break;

        case LGM_EXTMODEL_T01S:
                                m->Bfield = Lgm_B_T01S;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "T01S" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko 2001 Storm-time Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "References: (1) Tsyganenko, N. A., A model of the magnetosphere with a dawn-dusk asymmetry, 1, Mathematical structure, J. Geophys. Res., 107(A8), doi:10.1029/2001JA000219, 2002; (2) Tsyganenko, N. A., A model of the near magnetosphere with a dawn-dusk asymmetry, 2, Parameterization and fitting to observations, J. Geophys. Res., 107(A8), doi:10.1029/2001JA000220, 2002; (3) Tsyganenko, N. A., H. J. Singer, and J. C. Kasper, Storm-time distortion of the inner magnetosphere: How severe can it get? J. Geophys. Res., 108(A5), 1209, doi:10.1029/2002JA009808, 2003.");
                                strcpy( m->ExtMagModelStr4, "Comments: The T02 model optimized for storms. Only uses data from |X| < 15Re.");
                                break;

        case LGM_EXTMODEL_TS04:
                                m->Bfield = Lgm_B_TS04;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "TS04" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko Sitnov 2004 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A., and M. I. Sitnov (2005), Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms, J. Geophys. Res., 110, A03208, doi:10.1029/2004JA010798.");
                                strcpy( m->ExtMagModelStr4, "Comments: Also known as TS05.");
                                break;

        case LGM_EXTMODEL_TS07:
                                m->Bfield = Lgm_B_TS07;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "TS07" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko Sitnov 2007 Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "Reference: Tsyganenko, N. A., and M. I. Sitnov (2007), Magnetospheric configurations from a high-resolution data-based magnetic field model, J. Geophys. Res., 112, A06225, doi:10.1029/2007JA012260." );
                                strcpy( m->ExtMagModelStr4, "Comments: Coefficients are time-dependent. Coefficient files are required for this model to work properly.");
                                break;

        case LGM_EXTMODEL_TA16:
                                m->Bfield = Lgm_B_TA16;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "TA16" );
                                strcpy( m->ExtMagModelStr2, "Tsyganenko Andreeva 2016 RBF Magnetic Field Model" );
                                strcpy( m->ExtMagModelStr3, "References: Tsyganenko, N. A., and V. A. Andreeva (2016; doi: 10.1002/2016JA023217) and V. A. Andreeva and N. A. Tsyganenko (2016; doi: 10.1002/2015JA022242)" );
                                break;

        case LGM_EXTMODEL_OP77:
                                m->Bfield = Lgm_B_OP77;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "OP77" );
                                strcpy( m->ExtMagModelStr2, "Olson Pfitzer, 1977 Magnetic Field Model." );
                                strcpy( m->ExtMagModelStr3, "Reference: Olson, W. P., and K. A. Pfitzer (1977), Magnetospheric Magnetic Field Modeling, Ann. Sci. Rep. F44620-75-C-0033, Air Force Off.  of Sci. Res., Arlington, Va." );
                                strcpy( m->ExtMagModelStr4, "Comments: " );
                                break;

        case LGM_EXTMODEL_SCATTERED_DATA:
                                m->Bfield = Lgm_B_FromScatteredData;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;

        case LGM_EXTMODEL_SCATTERED_DATA2:
                                m->Bfield = Lgm_B_FromScatteredData2;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData2" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData2" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;

        case LGM_EXTMODEL_SCATTERED_DATA3:
                                m->Bfield = Lgm_B_FromScatteredData3;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData3" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData3" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;
        case LGM_EXTMODEL_SCATTERED_DATA4:
                                m->Bfield = Lgm_B_FromScatteredData4;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData4" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData4" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;
        case LGM_EXTMODEL_SCATTERED_DATA5:
                                m->Bfield = Lgm_B_FromScatteredData5;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData5" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData5" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;
        case LGM_EXTMODEL_SCATTERED_DATA6:
                                m->Bfield = Lgm_B_FromScatteredData6;
                                m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_RK5;
m->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
                                strcpy( m->ExtMagModelStr1, "ScatteredData6" );
                                strcpy( m->ExtMagModelStr2, "ScatteredData6" );
                                strcpy( m->ExtMagModelStr3, "Reference: Uses KDTree and nearest neighbor algorithm to interpolate from unstructured data clouds.");
                                strcpy( m->ExtMagModelStr4, "Comments: Any 3D collection of B-field data points can be used." );
                                break;


        default:
                                printf("Lgm_MagModelInfo_Set_MagModel(): No such B field model.\n");
                                exit(-1);
                                break;


    }

}

void Lgm_Get_IntMagModelStrings( char **Str1, char **Str2, char **Str3, char **Str4, Lgm_MagModelInfo *m ){
    *Str1 = m->IntMagModelStr1;
    *Str2 = m->IntMagModelStr2;
    *Str3 = m->IntMagModelStr3;
    *Str4 = m->IntMagModelStr4;
}

void Lgm_Get_ExtMagModelStrings( char **Str1, char **Str2, char **Str3, char **Str4, Lgm_MagModelInfo *m ){
    *Str1 = m->ExtMagModelStr1;
    *Str2 = m->ExtMagModelStr2;
    *Str3 = m->ExtMagModelStr3;
    *Str4 = m->ExtMagModelStr4;
}

void Lgm_Get_MagModel( int *InternalModel, int *ExternalModel, Lgm_MagModelInfo *m ){

    *InternalModel = m->InternalModel;
    *ExternalModel = m->ExternalModel;

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


void Lgm_Set_KdTree( Lgm_KdTree *KdTree, int k, double d2, Lgm_MagModelInfo *m ) {
    m->KdTree         = KdTree;
    m->KdTree_Alloced = TRUE;
    Lgm_Set_KdTree_kNN_k( m, k );
    Lgm_Set_KdTree_kNN_MaxDist2( m, d2 );
    return;
}
void Lgm_Set_KdTree_kNN_k( Lgm_MagModelInfo *m, int k ) {
    m->KdTree_kNN_k = k;
    return;
}
void Lgm_Set_KdTree_kNN_MaxDist2( Lgm_MagModelInfo *m, double MaxDist2 ) {
    m->KdTree_kNN_MaxDist2 = MaxDist2; // This maxdist is in physical units
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

void Lgm_Set_Lgm_B_TA16(Lgm_MagModelInfo *MagInfo) {
    MagInfo->Bfield = Lgm_B_TA16;
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

/*
 *  Set transformation matrix to go from GSM to PQB coords.
 *  
 *      B: z-like coordinate pointing in direction of local b.
 *      Q: y-like coordinate that is perp to both radial direction and b-hat
 *      P: x-like coordinate completes right handed system.
 */
void Lgm_Set_GSM_TO_PQB( Lgm_Vector *Position, Lgm_MagModelInfo *m ){

    Lgm_Vector P, Q, B, R;


    /*
     * Get B in GSM coords.
     */
    m->Bfield( Position, &B, m ); Lgm_NormalizeVector( &B ); // B in GSM coords

    /*
     *  The radial direction is just -Position
     */
    R.x = -Position->x; R.y = -Position->y; R.z = -Position->z; 

    /*
     * Compute  Q = R x B
     */
    Lgm_CrossProduct( &R, &B, &Q ); Lgm_NormalizeVector( &Q );   // Q in GSM coords
    
    /*
     * Compute  P = B x Q
     */
    Lgm_CrossProduct( &Q, &B, &P ); Lgm_NormalizeVector( &P );   // P in GSM coords

    m->Agsm_to_pqb[0][0] = P.x, m->Agsm_to_pqb[1][0] = P.y, m->Agsm_to_pqb[2][0] = P.z;
    m->Agsm_to_pqb[0][1] = Q.x, m->Agsm_to_pqb[1][1] = Q.y, m->Agsm_to_pqb[2][1] = Q.z;
    m->Agsm_to_pqb[0][2] = B.x, m->Agsm_to_pqb[1][2] = B.y, m->Agsm_to_pqb[2][2] = B.z;

    m->Apqb_to_gsm[0][0] = P.x, m->Apqb_to_gsm[1][0] = Q.x, m->Apqb_to_gsm[2][0] = B.x;
    m->Apqb_to_gsm[0][1] = P.y, m->Apqb_to_gsm[1][1] = Q.y, m->Apqb_to_gsm[2][1] = B.y;
    m->Apqb_to_gsm[0][2] = P.z, m->Apqb_to_gsm[1][2] = Q.z, m->Apqb_to_gsm[2][2] = B.z;

    m->Agsm_to_pqb_set = TRUE;

}

void Lgm_GSM_TO_PQB( Lgm_Vector *u_gsm, Lgm_Vector *u_pqb, Lgm_MagModelInfo *m ){

    if ( m->Agsm_to_pqb_set  == FALSE ){
        printf( "You need to set the transformations matrices first using Lgm_Set_GSM_TO_PQB()\n");
    } else {
        Lgm_MatTimesVec( m->Agsm_to_pqb, u_gsm, u_pqb );
    }

}


void Lgm_PQB_TO_GSM( Lgm_Vector *u_pqb, Lgm_Vector *u_gsm, Lgm_MagModelInfo *m ){

    if ( m->Agsm_to_pqb_set  == FALSE ){
        printf( "You need to set the transformations matrices first using Lgm_Set_GSM_TO_PQB()\n");
    } else {
        Lgm_MatTimesVec( m->Apqb_to_gsm, u_gsm, u_pqb );
    }

}


