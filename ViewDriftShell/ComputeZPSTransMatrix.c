#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Lgm/Lgm_CTrans.h>
#include <Lgm/Lgm_MagModelInfo.h>
#include "Rainbow2.h"

void ComputeZPSTransMatrix( Lgm_Vector *Rgeo, long int Date, double UT, double Agsm_to_ins[3][3], double Ains_to_gsm[3][3], double Qins_to_gsm[4], double Qgsm_to_ins[4] ) {

    double  Glon, Alpha, A, B, ca, sa, sum;
    Lgm_Vector  Xgeo, Ygeo, Zgeo, Xgsm, Ygsm, Zgsm;
    Lgm_Vector  Ngeo, Egeo, Ngsm, Egsm;
    Lgm_Vector  Xsc_geo,Ysc_geo, Zsc_geo;
    Lgm_Vector  Xsc_gsm,Ysc_gsm, Zsc_gsm;
    double  Agsm_to_sc[3][3], Asc_to_gsm[3][3], Asc_to_r3[3][3], Ar3_to_sc[3][3], Asc_to_dom[3][3], Adom_to_sc[3][3];  
    Lgm_Vector  X1, Y1, Z1;
    Lgm_Vector  X2, Y2, Z2;
    Lgm_Vector  X3, Y3, Z3;
    Lgm_Vector  X4, Y4, Z4;
    int     i, j, k;
    Lgm_CTrans  c;
    Lgm_Vector Xins, Yins, Zins, Xsc, Ysc, Zsc;
    Lgm_Vector Xins_new, Yins_new, Zins_new;
    double  Q1[4], Q2[4], Q3[4], QQ[4], QT1[4];



    Glon = atan2( Rgeo->y, Rgeo->x );



    /*
     *  Compute transformation matrix from dome coords to GSM
     *  Assume spacecraft is located at Rgeo.
     */
    
    
    /*
     * Nadir in GEO coords.
     */
    Ngeo.x = -Rgeo->x;
    Ngeo.y = -Rgeo->y;
    Ngeo.z = -Rgeo->z;
    Lgm_NormalizeVector( &Ngeo );
    

    
    /*
     *   East in GEO coords. Since, r = {cos(phi), sin(phi)}
     *   east direction will be dr/dphi which is;
     *
     *      east = {-sin(phi), cos(phi)}
     *  
     *   note that phi is just glon in GEO coords.
     *  
     */
    Egeo.x = -1.0*sin(Glon);
    Egeo.y = cos(Glon);
    Egeo.z = 0.0;


    /*
     * North in GEO coords.
     */
    Zgeo.x = 0.0;
    Zgeo.y = 0.0;
    Zgeo.z = 1.0;




    /*
     *  Define S/C coords as:
     *   X: East
     *   Y: South
     *   Z: Nadir
     */
    Xsc_geo = Egeo;
    Ysc_geo.x = -Zgeo.x; Ysc_geo.y = -Zgeo.y; Ysc_geo.z = -Zgeo.z;
    Zsc_geo = Ngeo;


    /*
     *  Convert to GSM
     */
    Lgm_Set_Coord_Transforms( Date, UT, &c );
    Lgm_Convert_Coords( &Xsc_geo, &Xsc_gsm, GEO_TO_GSM, &c );
    Lgm_Convert_Coords( &Ysc_geo, &Ysc_gsm, GEO_TO_GSM, &c );
    Lgm_Convert_Coords( &Zsc_geo, &Zsc_gsm, GEO_TO_GSM, &c );

Xsc_gsm = Xsc_geo;
Ysc_gsm = Ysc_geo;
Zsc_gsm = Zsc_geo;


    /*
     *  Construct Transformation Matrix from GSM to S/C
     */
    Agsm_to_sc[0][0] = Xsc_gsm.x, Agsm_to_sc[1][0] =  Xsc_gsm.y, Agsm_to_sc[2][0] = Xsc_gsm.z;
    Agsm_to_sc[0][1] = Ysc_gsm.x, Agsm_to_sc[1][1] =  Ysc_gsm.y, Agsm_to_sc[2][1] = Ysc_gsm.z;
    Agsm_to_sc[0][2] = Zsc_gsm.x, Agsm_to_sc[1][2] =  Zsc_gsm.y, Agsm_to_sc[2][2] = Zsc_gsm.z;

    Asc_to_gsm[0][0] = Xsc_gsm.x, Asc_to_gsm[1][0] =  Ysc_gsm.x, Asc_to_gsm[2][0] = Zsc_gsm.x;
    Asc_to_gsm[0][1] = Xsc_gsm.y, Asc_to_gsm[1][1] =  Ysc_gsm.y, Asc_to_gsm[2][1] = Zsc_gsm.y;
    Asc_to_gsm[0][2] = Xsc_gsm.z, Asc_to_gsm[1][2] =  Ysc_gsm.z, Asc_to_gsm[2][2] = Zsc_gsm.z;
printf("             ( %5g %5g %5g )\n", Asc_to_gsm[0][0],  Asc_to_gsm[1][0], Asc_to_gsm[2][0] );
printf("Asc_to_gsm = ( %5g %5g %5g )\n", Asc_to_gsm[0][1],  Asc_to_gsm[1][1], Asc_to_gsm[2][1] );
printf("             ( %5g %5g %5g )\n", Asc_to_gsm[0][2],  Asc_to_gsm[1][2], Asc_to_gsm[2][2] );




    /*
     *  Now construct trans matrix between S/C and INSTRUMENT
     *
     *  Lets define the instrument coord system to be as follows:
     *     Yins: parallel to instrument cylindrical symmetry axis.
     *     Zins: parallel to FOV of 
     */


    /*
     *  Here is the prescription to orient ZPS:
     *
     *      1. Start with Yins parallel to Ysc and Zins parallel to Zsc
     *      2. Rotate by -29 deg. about Zsc
     *      3. Then rotate by +45.5 deg. about Yins
     *      4. Then roate by +30.0 deg. about Zins
     * 
     * The easiest way to do this is to use Quaternions
     * 
     */
    
    // Step 1
    Xins = Xsc_gsm;
    Yins = Ysc_gsm;
    Zins = Zsc_gsm;

    // Step2
    //Lgm_AxisAngleToQuat( &Zins, -29.0, Q1 );
    Lgm_AxisAngleToQuat( &Zins, 0.0, Q1 );
    Lgm_QuatRotateVector( Q1, &Xins, &Xins_new ); Xins = Xins_new;
    Lgm_QuatRotateVector( Q1, &Yins, &Yins_new ); Yins = Yins_new;
    Lgm_QuatRotateVector( Q1, &Zins, &Zins_new ); Zins = Zins_new;

    // Step3
    //Lgm_AxisAngleToQuat( &Yins, 45.5, Q2 );
    Lgm_AxisAngleToQuat( &Yins, 0.0, Q2 );
    Lgm_QuatRotateVector( Q2, &Xins, &Xins_new ); Xins = Xins_new;
    Lgm_QuatRotateVector( Q2, &Yins, &Yins_new ); Yins = Yins_new;
    Lgm_QuatRotateVector( Q2, &Zins, &Zins_new ); Zins = Zins_new;
    
    // Step4
    //Lgm_AxisAngleToQuat( &Zins, 30.0, Q3 );
    Lgm_AxisAngleToQuat( &Zins, 0.0, Q3 );
    Lgm_QuatRotateVector( Q3, &Xins, &Xins_new ); Xins = Xins_new;
    Lgm_QuatRotateVector( Q3, &Yins, &Yins_new ); Yins = Yins_new;
    Lgm_QuatRotateVector( Q3, &Zins, &Zins_new ); Zins = Zins_new;

    /*
     *  Now, Xins, Yins, Zins will be oriented in their final positions.
     *  We can construct the sc -> ins trans matrix in one of two wyas now.
     *      1. combine all of the quats and then convert Q->matrix
     *  or  2. use the final unit vectors to manually construct the matri.
     *  Try both...
     */

    Lgm_QuatCombineQuats( Q1, Q2, QT1 );
    Lgm_QuatCombineQuats( QT1, Q3, QQ );
    for (i=0; i<4; i++) Qins_to_gsm[i] = QQ[i];

    Lgm_Quat_To_Matrix( QQ, Ains_to_gsm );

    for (i=0; i<3; ++i){ 
        for (j=0; j<3; ++j){ 
            Agsm_to_ins[i][j] = Ains_to_gsm[j][i];
        }
    }
    Lgm_MatrixToQuat( Agsm_to_ins, Qgsm_to_ins );

Lgm_MatrixToQuat( Agsm_to_sc, Qgsm_to_ins );
Lgm_MatrixToQuat( Asc_to_gsm, Qins_to_gsm );
}
