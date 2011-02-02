#ifndef LGM_FLUX_TO_PSD_H
#define LGM_FLUX_TO_PSD_H
#include "Lgm_DynamicMemory.h"
#include "Lgm_CTrans.h"
#include "Lgm_Vec.h"
static unsigned char Rainbow2_Red[] = { 0, 34, 35, 35, 36, 37, 38, 39, 39, 40, 40, 41,
41, 42, 42, 42, 42, 42, 42, 42, 41, 41, 41, 40, 40, 39, 38, 38, 37, 37, 36, 35,
35, 34, 33, 33, 32, 32, 31, 30, 30, 29, 28, 28, 27, 27, 26, 25, 25, 24, 24, 23,
23, 22, 22, 21, 21, 20, 20, 19, 19, 19, 18, 18, 17, 17, 17, 16, 16, 16, 15, 15,
15, 14, 14, 14, 13, 13, 13, 13, 12, 12, 12, 12, 11, 11, 11, 11, 11, 10, 10, 10,
10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 12, 17, 24, 29, 35,
42, 48, 53, 60, 67, 73, 79, 84, 91, 97, 103, 109, 115, 120, 126, 131, 136, 141,
147, 152, 156, 161, 165, 169, 173, 177, 181, 184, 188, 191, 194, 197, 200, 203,
206, 208, 211, 213, 215, 217, 220, 222, 224, 226, 228, 230, 232, 234, 236, 236,
236, 236, 236, 237, 237, 237, 237, 237, 238, 238, 238, 238, 238, 238, 239, 239,
239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 241, 241, 241, 241, 241, 241,
242, 242, 242, 242, 242, 242, 243, 243, 243, 243, 243, 243, 244, 244, 244, 244,
244, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 247, 247, 247, 247, 247,
248, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 251, 251, 251, 251, 252,
252, 252, 253, 253, 253, 253, 254, 254, 254, 255 };

static unsigned char Rainbow2_Grn[] = { 0, 44, 45, 47, 49, 51, 53, 55, 57, 58, 60, 62,
64, 66, 67, 69, 71, 72, 74, 75, 76, 78, 79, 80, 81, 82, 84, 85, 86, 87, 88, 90,
91, 92, 93, 95, 96, 98, 99, 101, 102, 104, 106, 107, 109, 111, 113, 115, 117,
119, 121, 124, 126, 129, 131, 134, 137, 140, 143, 146, 149, 152, 156, 160, 163,
167, 171, 176, 180, 185, 189, 194, 199, 205, 207, 208, 208, 209, 209, 210, 211,
211, 212, 212, 213, 213, 213, 214, 214, 215, 215, 216, 216, 217, 217, 217, 218,
218, 218, 219, 219, 219, 220, 220, 220, 221, 221, 221, 222, 222, 222, 223, 223,
223, 223, 224, 224, 224, 225, 225, 225, 225, 226, 226, 226, 227, 227, 227, 227,
228, 228, 228, 228, 229, 229, 229, 229, 230, 230, 230, 231, 231, 231, 231, 232,
232, 232, 232, 232, 233, 233, 233, 233, 234, 234, 234, 234, 234, 235, 235, 235,
235, 235, 236, 236, 234, 233, 231, 229, 228, 226, 225, 223, 221, 219, 217, 215,
213, 211, 209, 207, 205, 203, 201, 199, 198, 195, 194, 192, 190, 188, 186, 184,
183, 180, 178, 177, 175, 173, 171, 169, 167, 165, 163, 161, 159, 157, 156, 154,
151, 149, 147, 145, 143, 141, 139, 136, 134, 132, 129, 127, 125, 122, 119, 117,
114, 111, 108, 105, 102, 99, 95, 92, 89, 85, 82, 78, 74, 70, 67, 62, 58, 54,
50, 45, 41, 37, 33, 28, 24, 20, 16, 12, 7, 3, 0 };

static unsigned char Rainbow2_Blu[] = { 0, 89, 93, 96, 99, 102, 106, 109, 113, 116,
120, 123, 126, 130, 132, 136, 139, 141, 143, 146, 148, 150, 152, 154, 156, 158,
159, 161, 162, 164, 165, 167, 168, 169, 171, 172, 173, 174, 176, 177, 178, 179,
180, 181, 182, 184, 185, 186, 186, 188, 189, 189, 190, 191, 192, 193, 194, 195,
196, 196, 197, 198, 199, 199, 200, 201, 202, 202, 203, 204, 205, 205, 206, 207,
205, 200, 197, 192, 187, 183, 178, 173, 168, 162, 157, 153, 147, 141, 135, 130,
124, 118, 112, 106, 100, 95, 89, 83, 77, 71, 65, 60, 54, 48, 43, 38, 32, 26,
21, 14, 9, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


#define LGM_Ee0     0.510998910  // Electron rest energy in MeV

typedef unsigned char byte; // needed for gif writer

typedef struct Lgm_FluxToPsd {


    /*
     * Params for Differential Flux versus Energy and Pitch Angle.
     * And for PSD versus Energy and Pitch Angle.
     */
    int          nE;                //!< Number of energy bins in Flux array.

    double       *E;                //!< Array of energy values in Flux array.

    int          nA;                //!< Number of pitch angle bins in Flux array.

    double       *A;                //!< Array of pitch angle values in Flux array.

    double       **FLUX_EA;         //!< Array of differential flux versus Energy and PitchAngle, Flux[E][A].

    double       **PSD_EA;          //!< Array of PSD versus Energy and PitchAngle, PSD[E][A].

    int          Alloced1;          //!< If true, the arrays are alloced.



    /*
     * Intermediate quantities needed.
     */
    Lgm_DateTime DateTime;         //!< Date/Time of measurment.

    Lgm_Vector   Position;         //!< Position of measurment.

    double       *AofK;            //!< Array of Alpha values that are implied by the k values. Size is nK.

    double       **EofMu;          //!< Array of Energy values that are implied by the Mu, Alpha and B values. Size is nMu.

    double       B;                //!< Magnetic field strength.


    /*
     * Params for PSD versus Mu and K
     */
    int          nMu;               //!< Number of Mu bins in PSD array.

    double       *Mu;               //!< Array of Mu values in PSD array.

    int          nK;                //!< Number of K bins in PSD array.

    double       *K;                //!< Array of K values in PSD array.


    double       **PSD_MK;          //!< Array of PSD versus Mu and K,  PSD[Mu][K].

    int          Alloced2;          //!< If true, the arrays are alloced.

    




    /*
     * Other things..
     */
    int          DumpDiagnostics;    //!< If true, some diagnostics (images, etc) may get dumped out.



} Lgm_FluxToPsd;



Lgm_FluxToPsd *Lgm_CreateFluxToPsd( int DumpDiagnostics );
void           Lgm_FreeFluxToPsd( Lgm_FluxToPsd *f );
void           Lgm_FluxToPsd_SetFlux( double **J, double *E, int nE, double *A, int nA, Lgm_FluxToPsd *f );
void           Lgm_FluxToPsd_SetDateTimeAndPos( Lgm_DateTime *d, Lgm_Vector *u, Lgm_FluxToPsd *f );
void           Lgm_FluxPsd_GetPsdAtConstMusAndKs( double *Mu, int nMu, double *K, int nK, Lgm_FluxToPsd *p );
double         Lgm_FluxPsd_GetPsdAtEandAlpha( double E, double a, Lgm_FluxToPsd *f );






void   DumpGif( char *Filename, int W, int H, double **Image );
double Lgm_Energy_to_Mu( double E, double a, double B );
double Lgm_Mu_to_Energy( double Mu, double a, double B );
double Lgm_p2c2( double Ek, double E0 );
double Lgm_v2overc2( double Ek, double E0 );
double Lgm_gamma( double Ek, double E0 );
double Lgm_PsdToDiffFlux( double f, double p2c2 );
double Lgm_DiffFluxToPsd( double j, double p2c2 );



#endif
