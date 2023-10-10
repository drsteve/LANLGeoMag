#include <Lgm_MagEphemInfo.h>
#include <stdio.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
#include <Lgm_ElapsedTime.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define ELECTRON 0
#define PROTON   1

#define GCA_ONLY          0
#define FULL_ORBIT_ONLY   1
#define ADAPTIVE          2

#define GCA         0
#define FULL_ORBIT  1

typedef struct Lgm_ParticleInfo {

    Lgm_MagModelInfo    *mInfo;


    /*
     * Particle properties
     */
    double              q;          //<! particle charge
    double              m;          //<! particle mass
    double              mu;         //<! relativistic magnetic moment
    double              KE;         //<! initial kinetic energy (MeV)
    double              AlphaEq;    //<! initial equatorial pitch angle (Degrees)
    double              Alpha;      //<! initial pitch angle (Degrees)
    double              B0;         //<! initial |B| (A)nT
    double              c, c2;      //<! speed of light, and its square

    double              mc, mc2;    //<! mc amd mc^2
    double              gamma;      //<! relativistic factor


    int                 CurrentIntegrator;

    /*
     * X and v_par (and p_par)
     */
    Lgm_Vector          X;          //<! Guiding Center position
    double              v_par;      //<! parallel velocity
    double              p_par;      //<! parallel momentum

    /*
     *  Boris velocities
     */
    Lgm_Vector u;                   //<! staggered velocity used in Boris scheme


    /*
     * x and v
     */
    Lgm_Vector          x;          //<! Full Orbit position
    Lgm_Vector          v;          //<! Full Orbit velocity



    /*
     * ODE integrators
     */
    gsl_odeiv2_system Sys_Full;
    gsl_odeiv2_driver *Driver_Full;

    gsl_odeiv2_system Sys_GCA;
    gsl_odeiv2_driver *Driver_GCA;










    /*
     * Derived results
     */
    Lgm_Vector          Es;         //<! Estar
    Lgm_Vector          Bs;         //<! Bstar
    Lgm_Vector          Es_cross_b; //<! Es x b
    double              Bs_par;     //<! parallel comp of Bs
    Lgm_Vector          E;          //<! E
    Lgm_Vector          B;          //<! B
    double              Bmag;       //<! |B|
    Lgm_Vector          b;          //<! B-hat
    Lgm_Vector          dbdt;       //<! dbdt
    Lgm_Vector          Grad_B;     //<! Gradient of B
    Lgm_Vector          Grad_b;     //<! Gradient of b
    Lgm_Vector          Curl_B;     //<! Curl of B
    Lgm_Vector          Curl_b;     //<! Curl of b
    Lgm_Vector          Curl_E;     //<! Curl of E
    Lgm_Vector          Grad_gamma; //<! Gradient of Relativistic factor.

    Lgm_Vector          X_dot;      //<! velocity of guiding center position
//    double              v_par_dot;  //<! rate of change of v_par
    double              p_par_dot;  //<! rate of change of v_par

} Lgm_ParticleInfo;


/*
 * Function prototypes
 */
int Lgm_TraceParticle_GCA( Lgm_Vector *ParticlePosition, long int *nParticlePosition  );
int Lgm_TraceParticle_Boris( Lgm_Vector *ParticlePosition, long int *nParticlePosition  );
int Lgm_TraceParticle_Full( Lgm_Vector *ParticlePosition, long int *nParticlePosition,
Lgm_Vector *ParticlePosition2, long int *nParticlePosition2  );

int InitParticleIntegrators( Lgm_ParticleInfo *p );




