/**
 *  Returns a pointer to a dynamically allocated Lgm_PsdToFlux structure.
 *  User must destroy this with Lgm_P2F_FreePsdToFlux() when done.
 *
 *      \param[in]      DumpDiagnostics  Boolean flag to turn on/off dumping of diagnostics.
 *      \return         A pointer to an allocated and initialized Lgm_PsdToFlux stucture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
Lgm_FluxToPsd *Lgm_P2F_CreatePsdToFlux( int DumpDiagnostics ) {

    Lgm_FluxToPsd *p;

    /*
     * Allocate memory for a Lgm_PsdToFlux structure.
     */
    p = (Lgm_PsdToFlux *) calloc( 1, sizeof(*p) );

    /*
     * Set DumpDiagnostics flag to what we got here. This can be changed later as well.
     */
    p->DumpDiagnostics = DumpDiagnostics;

    p->Extrapolate = TRUE;

    p->Alloced1 = FALSE;
    p->Alloced2 = FALSE;


    return p;

}

/**
 * Destroy a dynamically allocated Lgm_PsdToFlux structure. (E.g. one that was
 * created by Lgm_P2F_CreatePsdToFlux().)
 *
 *      \param          p  Pointer to the allocated Lgm_PsdToFlux structure that you want to destroy.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_FreePsdToFlux( Lgm_PsdToFlux *p ) {

    if ( p->Alloced1 ) {
        LGM_ARRAY_1D_FREE( p->E );
        LGM_ARRAY_1D_FREE( p->A );
        LGM_ARRAY_2D_FREE( p->FLUX_EA );
        LGM_ARRAY_2D_FREE( p->PSD_EA );
    }

    if ( p->Alloced2 ) {
        LGM_ARRAY_1D_FREE( p->Mu );
        LGM_ARRAY_1D_FREE( p->K );
        LGM_ARRAY_1D_FREE( p->AofK );
        LGM_ARRAY_2D_FREE( p->EofMu );
        LGM_ARRAY_2D_FREE( p->PSD_MK );
    }

    free( p );

    return;
}


/**
 *  Set Date/Time and position in the Lgm_PsdToFlux structure.
 *      
 *     
 *      \param[in]      d   Date/Time of measurement.
 *      \param[in]      u   Position of measurment (in GSM).
 *      \param[in,out]  p   Lgm_PsdToFlux sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_SetDateTimeAndPos( Lgm_DateTime *d, Lgm_Vector *u, Lgm_PsdToFlux *p ) {

    p->DateTime = *d;
    p->Position = *u;

}


/**
 *     Adds (to a Lgm_PsdToFlux structure) the user-supplied arrays containing PSD[Mu][K],  Mu[], K[]
 *
 *      \param[in]      P                 2D array containing the Phase Space Density as a function of Mu and K.
 *      \param[in]      Mu                1D array containing the energy values implied by the first index of Flux[][] array.
 *      \param[in]      nMu               number of energies.
 *      \param[in]      K                 1D array containing the pitch angles values implied by the second index of Flux[][] array.
 *      \param[in]      nK                number of pitch angles.
 *      \param[in,out]  p                 Lgm_FluxToPsd sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_SetPsd( double **P, double *Mu, int nMu, double *K, int nK, Lgm_PsdToFlux *p ) {

    
    int     i, j;


    /*
     * If arrays are already alloc'd, free them first
     */
    if ( p->Alloced1 ) {
        LGM_ARRAY_1D_FREE( p->Mu );
        LGM_ARRAY_1D_FREE( p->K );
        LGM_ARRAY_2D_FREE( p->PSD_MK );
    }


    /*
     * Add Psd array to p structure. Alloc arrays appropriately.
     */
    p->nMu = nMu; 
    p->nK  = nK; 
    LGM_ARRAY_1D( p->Mu, p->nMu, double );
    LGM_ARRAY_1D( p->K, p->nK, double );
    LGM_ARRAY_2D( p->PSD_MK, p->nMu, p->nK, double );
    for (i=0; i<p->nMu; i++) p->Mu[i] = Mu[i];
    for (i=0; i<p->nK; i++)  p->K[i]  = K[i];
    for (i=0; i < p->nMu; i++) {
        for (j=0; j < p->nK; j++) {
            p->PSD_MK[i][j] = P[i][j]; // PSD_MK is "PSD versus Mu and K".
        }
    }
    if ( p->DumpDiagnostics ) {
        DumpGif( "Lgm_Psd_ToFlux_SetPsd_PSD_MK", f->nK, MK>nMu, MK>PSD_MK );
    }


    f->Alloced1 = TRUE;
   
    return;

}




/**
 *  Computes Flux at user-supplied constant values of E
 *  and \f$\alpga\f$.
 *
 *  This routine ( Lgm_P2F_GetFluxAtConstEsAndAs() ) must operate on a
 *  pre-initialized Lgm_PsdToFlux structure.  The routine Lgm_P2F_SetPsd()
 *  is used to add PSD data/info to an Lgm_PsdToFlux structure.
 *
 *  We want Flux at constant E and \f$\alpha \f$ (i.e. we want \J$f( E, \alpha
 *  )\f$).  
 * 
 *  To perform the calculation we note that \f$f( \mu, K )\f$ is the same as
 *  \f$f( \mu(E, \alpha), K(\alpha) )\f$. Thus, for a given E and \f$\alpha\f$,
 *  we can figure out what \f$\mu\f$ and K  they correspond to and then we can
 *  just use the \f$f(\mu, K)\f$ array to compute the desired f values. Then
 *  finally convert f to J.
 *  The steps are;
 *
 *      - For each \f$\alpha\f$ compute \f$K(\alpha)\f$. This is done with the routine
 *        Lgm_KofAlpha().
 *
 *      - Then we compute \f$\mu\f$ from E and \f$\alpha\f$.
 *
 *      - Then we look up \f$f(\mu, K)\f$ from the array (interp or fit or
 *        whatever).
 *
 *      \param[in]      nE          Number of E values
 *      \param[in]      E           1-D array of E values
 *      \param[in]      nA          Number of Alpha values
 *      \param[in]      A           1-D array of Alpha values
 *      \param[in]      Extrapolate Turns on/off extrapolation capability
 *      \param[in,out]  p           A properly pre-initialized Lgm_PsdToFlux structure.
 *
 *      \author     Mike Henderson
 *      \date       2011
 *      \warning    Still working on this code. It is not finished.
 * 
 */
void Lgm_P2F_GetFluxAtConstEsAndAs( double *E, int nE, double *A, int nA, Lgm_PsdToFlux *p ) {

int                 k, m, DoIt;
double              AlphaEq, SinA;
Lgm_MagModelInfo    *mInfo, *mInfo2;

    /*
     * Init mInfo
This is no good! How does user define mag model etc...?
I think there needs to be a Lgm_MagModelInfo struct in f
Then add a routine to set stuff up in there. Or just use the ones we have already....
For now we will just go with the defaults.
     */
    mInfo = Lgm_InitMagInfo();
//mInfo->Bfield = Lgm_B_edip;

    /*
     * If arrays are already alloc'd, free them first
     */
    if ( p->Alloced2 ) {
        LGM_ARRAY_1D_FREE( p->E );
        LGM_ARRAY_1D_FREE( p->A );
        LGM_ARRAY_1D_FREE( p->KofA );
        LGM_ARRAY_2D_FREE( p->MuofE );
        LGM_ARRAY_2D_FREE( p->PSD_EA );
        LGM_ARRAY_2D_FREE( p->FLUX_EA );
    }
    
    /*
     * Alloc arrays
     */
    p->nE = nE; 
    p->nA = nA; 
    LGM_ARRAY_1D( p->E,     p->nE, double );
    LGM_ARRAY_1D( p->A,     p->nA,  double );
    LGM_ARRAY_1D( p->AofK,  p->nA,  double );
    LGM_ARRAY_2D( p->MuofE, p->nE,  p->nA,  double );


    /*
     * Copy A's (given in the arguments) into p structure.
     * Transform the A's into K's using Lgm_KofAlpha().
     * Save the results in the p structure.
     */
    Lgm_Setup_AlphaOfK( &(p->DateTime), &(p->Position), mInfo );
    p->B = mInfo->Blocal;
    {
        #pragma omp parallel private(mInfo2,SinAlphaEq,AlphaEq)
        #pragma omp for schedule(dynamic, 1)
        for ( k=0; k<nA; k++ ){

            mInfo2 = Lgm_CopyMagInfo( mInfo );  // make a private (per-thread) copy of mInfo

            p->A[k]    = A[k]; // A is local Pitch Angle
            SinAlphaEq = sqrt( mInfo2->Bmin/mInfo2->Blocal ) * sin( RadPerDeg*p->A[k] );
            AlphaEq    = DegPerRad*asin( SinAlphaEq );
            p->KofA[k] = Lgm_KofAlpha( AlphaEq, mInfo2 );

            Lgm_FreeMagInfo( mInfo2 ); // free mInfo2

        }
    }
    Lgm_TearDown_AlphaOfK( mInfo );


    /*
     * Copy E's (given in the arguments) into p structure.
     * Transform the E's into Mu's.
     * Save the results in the p structure.
     * Note that since this conversion involves E and Alpha, the result is 2D.
assumes electrons -- generalize this...
     */
    for ( m=0; m<nE; m++ ){
        p->E[m] = E[m];
        for ( k=0; k<nA; k++ ){
            p->MuofE[m][k] = Lgm_Ek_to_Mu( p->E[m], p->A[k], p->B, LGM_Ee0 );
            //printf("f->Mu[%d], f->K[%d], f->AofK[%d], f->B, f->EofMu[%d][%d] = %g %g %g %g %g\n", m, k, k, m, k, f->Mu[m], f->K[k], f->AofK[k], f->B, f->EofMu[m][k]);
        }
    }


    /*
     * Now, from the PSD[Mu][K] array, get PSD at the Mu's and K's we just computed.
     * The result will be the same as PSD at the given E's and A's
     */
    LGM_ARRAY_2D( p->PSD_EA,  p->nE,  p->nA,  double );
    LGM_ARRAY_2D( p->FLUX_EA, p->nE,  p->nA,  double );
    for ( m=0; m<nE; m++ ){
        for ( k=0; k<nA; k++ ){
            DoIt = FALSE;

            if ( p->Extrapolate > 2 ){ // extrapolate above and below

                DoIt = TRUE;

            } else if ( p->Extrapolate == 2) {

                if (p->MuofE[m][k] >= p->Mu[0] ) DoIt = TRUE; // extrapolate above

            } else if ( p->Extrapolate == 1) {

                if (p->MuofE[m][k] <= p->Mu[p->nMu-1] ) DoIt = TRUE; // extrapolate below

            } else if ( p->Extrapolate == 0) {

                if ((p->MuofE[m][k] >= p->Mu[0])&&(p->MuofE[m][k] <= p->Mu[p->nMu-1])) DoIt = TRUE; // interp only

            }

            if (DoIt) {
                p->PSD_EA[m][k]  = Lgm_P2F_GetPsdAtMuAndK( p->MuofE[m][k], p->KofA[k], p->A[k], p );
                // Now do conversion from units of PSD to Flux
                p2c2 = Lgm_p2c2( p->E[m], LGM_Ee0 );
                p->FLUX_EA[m][k] = Lgm_PsdToDiffFlux( p->PSD_EA[m][k], p2c2 );
            } else {
                p->PSD_EA[m][k] = 0.0;
            }

        }
    }

    if ( p->DumpDiagnostics ) {
        DumpGif( "Lgm_PsdToFlux_SetPsd_PSD_EA", p->nA, p->nE, p->PSD_EA );
    }


// FIX this... The create and destroy of this should not be in here...
    Lgm_FreeMagInfo( mInfo );

    p->Alloced2 = TRUE;

    return;

}



/**
 * The p structure should have an initialized PSD[Mu][K] array in it (i.e. as
 * added by Lgm_P2F_SetPsd()).  This routine computes psd from this array given
 * arbitrary values of Mu and K. (Alpha is also needed here, because we need to
 * convert some Mu's to to Energie's long the way).
 */
double  Lgm_P2F_GetPsdAtMuAndK( double Mu, double K, double A, Lgm_PsdToFlux *p ) {

    int         j, i, i0, i1;
    double      K0, K1, y0, y1, slp, psd;
    _FitData    *FitData;

    // if K < 0, we should return fill value.
    if ( K < 0.0 ) return(-9e99);

    FitData = (_FitData *) calloc( 1, sizeof( _FitData ) );


    /*
     * Interpolate on K first to get a 1D array of f(mu).
     */
    if ( K < p->K[0] ) {
        return(-9e99);
        i0 = 0; i1 = 1;
    } else if ( K > p->K[p->nK - 1] ) {
        return(-9e99);
        i0 = p->nK - 2; i1 = p->nK - 1;
    } else {
        for (i=1; i<p->nK; i++) {
            if ( K < p->K[i] ) {
                i0 = i-1; i1 = i;
                break;
            }
        }
    }
    //printf("i0, i1 = %d %d\n", i0, i1);


    // interpolate K
    FitData->n = p->nMu;
    LGM_ARRAY_1D( FitData->E, FitData->n, double );
    LGM_ARRAY_1D( FitData->g, FitData->n, double );
    for (j=0; j<f->nE; ++j){
        K0   = p->K[i0];
        K1   = p->K[i1];
        y0   = p->PSD_MK[j][i0];
        y1   = p->PSD_MK[j][i1];
        slp  = (y1-y0)/(K1-K0);
        FitData->g[j] = slp*(K-K1) + y1;
        if (FitData->g[j] < 0.0) FitData->g[j] = y0; // probably shouldnt be in the fit at all!?
        FitData->E[j] = Lgm_Mu_to_Ek( p->Mu[j], A, p->B, LGM_Ee0 );
//        printf("i0, i1 = %d %d   a0, a1 = %g %g   y0, y1, slp = %g %g %g    a = %g, FitData->g[%d] = %g\n", i0, i1, a0, a1, y0, y1, slp, a, j, FitData->g[j]);
    }


        
    // interpolate/fit E
    // for now just do a linear interp.
    // no lets try a fit...
    double  in[10], out[7], x[5];
    in[0] = 1e-8;
    in[1] = in[2] = 1e-9; //Info->Praxis_Tolerance;
    in[5] = 30000.0; //(double)Info->Praxis_Max_Function_Evals;
    in[6] = 10.0; //Info->Praxis_Maximum_Step_Size;
    in[7] = 10.0; //Info->Praxis_Bad_Scale_Paramater;
    in[8] = 4.0; //(double)Info->Praxis_Max_Its_Without_Improvement;
    in[9] = 1.0; //(double)Info->Praxis_Ill_Conditioned_Problem;
    x[0] = 0.0;
    x[1] = -1.0;
    x[2] = 25.0;
    x[3] = -2.0;
    x[4] = 200.0;
    praxis( 4, x, (void *)FitData, Cost, in, out);
/*
printf("out[0] = %g\n", out[0]);
printf("out[1] = %g\n", out[1]);
printf("out[2] = %g\n", out[2]);
printf("out[3] = %g\n", out[3]);
printf("out[4] = %g\n", out[4]);
printf("out[5] = %g\n", out[5]);
printf("out[6] = %g\n", out[6]);
printf("x[1] = %g   x[2] = %g   Cost = %g\n", x[1], x[2], out[6]);

FILE *fp;
printf("E = %g\n", E);
fp = fopen("data.txt", "w");
for (j=0; j<f->nE; ++j){
fprintf(fp, "%g %g\n", f->E[j], FitData->g[j]);
}
fclose(fp);
    
exit(0);
*/
    

//x[2] = 200.0;
//printf("x - %g %g %g %g\n", x[1], x[2], x[3], x[4]);
    E = Lgm_Mu_to_Ek( Mu, A, p->B, LGM_Ee0 );
    psd = Model( x,  E );
    //psd = (double)a;

//printf("E, a = %g %g  x = %g %g psd = %g\n", E, a, x[1], x[2], psd);
    
    

    LGM_ARRAY_1D_FREE( FitData->E );
    LGM_ARRAY_1D_FREE( FitData->g );
    free( FitData );


    return( psd );

}

