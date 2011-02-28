#ifndef PSD_ASSIM_H
#define PSD_ASSIM_H

#define ELECTRON    0
#define PROTON      1

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

typedef struct {

    /*
     * A consistent label for the B-field model used.
     */
    char    FieldModelStr[80];

    /*
     * local magnitude of B
     */
    double  Blocal;
    double  Blocal_mod;
    double  Blocal_obs;

    /*
     * Keep track of dipole moments (important for cross-comparing L* calcs)
     */
    double  Mused; // value of M that we used to compute L*
    double  Mref;  // A standard M value...
    double  Migrf; // true (time-dep.) value of M

    /*
     *  Define the K values that we want to use.
     *  Also need vars for corresponding quantities.
     */
    int     nK;
    double  K[5];
    double  Alpha[5];
    double  Bm[5];
    double  I[5];
    double  Lstar[5];


    /*
     *  Define the mu values that we want to use.
     */
    int     nMu;
    double  Mu[5];


    /*
     *  Define an array to hold the Mu versus Alpha PSD Image
     *  that we need to do interpolation in.
     */
    double  **PSD_vs_Mu_and_Alpha;
    double  *PSD_Mus;
    double  *PSD_Alphas;


    /*
     * Define an nMu x nK array to hold PSD vals
     */
    double  PSD[5][5];


} _PsdAssim;


/*
 * Protos
 */




/*
 * Global defs
 */
#ifdef MAIN
/*
 *  Make PsdAssim Global. We could get around this,
 *  but it adds complexity. Also be careful with threads.
 */
_PsdAssim   *PsdAssim;
#else
extern _PsdAssim   *PsdAssim;
#endif


#endif
