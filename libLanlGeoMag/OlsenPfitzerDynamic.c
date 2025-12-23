#include "Lgm/Lgm_MagModelInfo.h"

/*
 *   Olson-Pfitzer 1988 Dynamic Model.
 *   Converted to thread-safe C by M. Henderson 2013.
 *
 */
void Lgm_OP88_BDYN( double DEN, double VEL, double DST, double X, double Y, double Z, double *BX, double *BY, double *BZ )  {

    double  XX[4], B[4], SOFFD, SRING, STAIL, STDOFF, RINGST;

    SOFFD = Lgm_OP88_STDOFF( VEL, DEN );
    SRING = Lgm_OP88_RINGST( SOFFD, DST );
    STAIL = 1.0;

    // GET THE EXTERNAL CONTRIBUTION
    XX[1] = X; XX[2] = Y; XX[3] = Z;
    Lgm_OP88_BDYNAM( XX, B, SOFFD, SRING, STAIL );
    *BX = B[1]; *BY = B[2]; *BZ = B[3];

    return;

}


/*
 *  VERSION 5/13/88
 *  DEVELOPED MCDONNELL DOUGLAS
 *  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
 *
 *     PURPOSE
 *        CALCULATE THE TOTAL EXTERNAL MAGNETIC FIELD DURING DISTURBED TIMES.
 *
 *     METHOD
 *        CALLS THE EXTERNAL QUIET TIME SUBROUTINES AND COMBINES THEM
 *        ACCORDING TO THE DYNAMIC SCALING ALGORITHMS.
 *
 *     INPUT -- ARGUMENT LIST
 *        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD
 *               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY
 *               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
 *               Z IS ALONG THE EARTH'S NORTH DIPOLE AXIS. X IS PERPENDICULAR
 *               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
 *               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR
 *               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM.
 *
 *        SOFFD  THE STANDOFF DISTANCE OF THE MAGNETOPAUSE. THE QUIET STANDOFF
 *               DISTANCE IS 10.5 EARTH RADII. ACCEPTABLE VALUES RANGE BETWEEN
 *               6 AND 11. THIS VALUE IS USED TO CALCULATE THE STRENGTH OF
 *               THE MAGNETOPAUSE CURRENTS AND TO SCALE THE SIZE OF THE
 *               MAGNETOPAUSE. THIS VALUE ALSO SCALES THE SIZE OF THE TAIL
 *               CURRENT SYSTEM. THE RING SYSTEM IS NOT SCALED, SINCE ITS
 *               SOURCE IS PRIMARILY AT RADIAL DISTANCES.
 *
 *        SRING  RELATIVE STRENGTH OF THE RING CURRENT. A VALUE OF ONE
 *               UTILIZES THE NOMINAL QUIET RING VALUES BUILT INTO THE BASIC
 *               MODEL. THIS BASIC MODEL HAS A MAXIMUM RING DEPRESSION OF
 *               40 NT AT L=4 RE. IF SRING IS SET TO 2 THE RING DEPRESSION
 *               WOULD BE 80 NT.
 *
 *        STAIL  A TAIL CURRENT STRENGTH MULTIPLIER. WHEN STAIL IS EQUALTO
 *               1.0 THEN THE TAIL SCALES WITH THE STRENGTH OF THE MAGNETOPAUSE
 *               CURRENTS. TO WEAKEN THE TAIL FROM THIS VALUE USE VALUES
 *               LESS THAN 1.0, TO STRENGTHEN USE VALUES GREATER THAN 1.0.
 *
 *     OUTPUT -- ARGUMENT LIST
 *        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE
 *               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.
 *               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
 *               THE UNITS ARE GAUSS.
 */
void Lgm_OP88_BDYNAM( double *XX, double *BB, double SOFFD, double SRING, double STAIL ) {

    double  SCL, STRMAG, STRRIN, STRTAI, XXX[4], XXXX[4], BM[4], BR[4], BT[4];
    int     IDX;

    // CALCULATE MAGNETOPAUSE SCALE FACTOR.
    SCL = 10.5/SOFFD;

    // CALCULATE STRENGTH OF MAGNETOPAUSE CURRENTS.
    STRMAG = SCL*SCL*SCL;

    // SET STRENGTH OF RING AND TAIL.
    STRRIN = SRING;
    STRTAI = STAIL*STRMAG;

    // CALCULATE SCALED DISTANCES
    for ( IDX=1; IDX<=3; IDX++ ) {
        XXX[IDX]  = XX[IDX]*SCL;
        XXXX[IDX] = XX[IDX];
    }

    // CALL THE QUIET TIME SUBROUTINE.
    Lgm_OP88_BFMAGP( XXX,  BM );
    Lgm_OP88_BFRING( XXXX, BR );
    Lgm_OP88_BFTAIL( XXX,  BT );

    /*
     *   COMBINE THE COMPONENTS OF THE MAGNETIC FIELD ACCORDING TO THEIR
     *   RELATIVE STRENGTHS.
     */
    for ( IDX=1; IDX<=3; IDX++ ) BB[IDX] = STRMAG*BM[IDX] + STRRIN*BR[IDX] + STRTAI*BT[IDX];

    return;

}



/*
 *  VERSION 5/13/88
 *  DEVELOPED MCDONNELL DOUGLAS
 *  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
 *
 *     INPUT -- ARGUMENT LIST
 *        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD
 *               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY
 *               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
 *               Z IS ALONG THE EARTH'S NORTH DIPOLE AXIS. X IS PERPENDICULAR
 *               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
 *               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR
 *               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM.
 *
 *     OUTPUT -- ARGUMENT LIST
 *        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE
 *               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.
 *               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
 *               THE UNITS ARE NANOTESLA.
 *
 * THIS SUBROUTINE CONTAINS NEW, REFITTED COEFFICIENTS FOR COMPUTING ALL THE
 * B-MAGNETOPAUSE COMPONENTS.
 * THE FORM OF THE EXPANSION IS GIVEN IN THE NEXT FEW STATEMENTS.
 *
 * BMPX=(1 + X + X**2 + X**3 + X**4 + 1/(10-(X-2)**2))*(1 + Y**2 + Y**4)*
 *      (Z + Z**3 + Z**5).
 * BMPY=(1 + X + X**2 + X**3 + X**4)*(Y + Y**3 + Y**5)*(Z + Z**3 + Z**5).
 * BMPZ=(1 + X + X**2 + X**3 + X**4 + 1/(15-X) + 1/((30-X)**2))*
 *      (1 + Y**2 + Y**4)*(1 + Z**2 + Z**4).
 *
 * COEFFICIENTS COMPUTED FROM COMBINED OLSON DATA-SETS BOUNDARY.DAT AND
 * BOUND.DAT, FOR Y>0 AND Z>0 ONLY, I.E. A TOTAL OF 1009 DATA POINTS.
 * 2 EXTRA EXTRAPOLATED POINTS ADDED FOR Z-COEFF, TO IMPROVE FIT,
 * NAMELY X=10,Y=Z=0,BZ=29 AND X=11,Y=Z=0,BZ=30.25.
 *
 * ***PRELIMINARY ROUTINE
 * ***VALID TO APPROXIMATELY -60 RE
 * ***MAGNETIC FIELD FROM MAGNETOPAUSE CURRENTS ONLY
 */
void Lgm_OP88_BFMAGP( double *XX, double *BB ) {

    double  XN, YN, ZN, X[6], Y[6], Z[6], g, g2, FX1, XF1, XF2;
    int     IDX;


    static double A[] = {   -9e99,  0.113275039e+01,  0.354408138e-01, -0.152252289e-02,
                                   -0.683306571e-04, -0.642841428e-06, -0.121504674e-01,
                                   -0.839622808e-03, -0.167520029e-04, -0.385962942e-07,
                                    0.107674747e-08,  0.558984066e-04,  0.551508083e-05,
                                    0.206288036e-06,  0.335316730e-08,  0.198413126e-10,
                                   -0.545824692e-02, -0.264107861e-03,  0.143533146e-05,
                                    0.195177861e-06,  0.207546358e-08,  0.211199178e-03,
                                    0.220245929e-04,  0.860991804e-06,  0.145349395e-07,
                                    0.886173426e-10, -0.949615014e-06, -0.110830563e-06,
                                   -0.477998707e-08, -0.873645670e-10, -0.569051859e-12,
                                    0.271760982e-04,  0.266707661e-05,  0.994617153e-07,
                                    0.167023062e-08,  0.104617062e-10, -0.989193381e-06,
                                   -0.113236254e-06, -0.482686247e-08, -0.880319914e-10,
                                   -0.575385009e-12,  0.487020380e-08,  0.586310778e-09,
                                    0.260182431e-10,  0.488435735e-12,  0.326678627e-14,
                                    0.193470073e+01,  0.402453184e+00, -0.193471275e-02,
                                    0.682263300e-01, -0.576195028e-02,  0.237557251e-04,
                                   -0.529665092e-03,  0.255710365e-04, -0.120115033e-06 };

    static double B[] = {  -9e99,  -0.519952811e-01, -0.230140495e-02,  0.146173188e-03,
                                    0.809832090e-05,  0.888401672e-07, -0.370911323e-03,
                                   -0.101231737e-03, -0.742647399e-05, -0.196170248e-06,
                                   -0.165503899e-08,  0.150949325e-05,  0.308240260e-06,
                                    0.195390104e-07,  0.472441419e-09,  0.375989214e-11,
                                   -0.422217818e-04, -0.621468353e-04, -0.620102765e-05,
                                   -0.189322407e-06, -0.172039538e-08,  0.445292017e-05,
                                    0.118324999e-05,  0.855768008e-07,  0.223059815e-08,
                                    0.183677951e-10, -0.550030643e-08, -0.150351465e-08,
                                   -0.107031245e-09, -0.268793755e-11, -0.205845354e-13,
                                    0.813519478e-06,  0.279971147e-06,  0.227601529e-07,
                                    0.643000209e-09,  0.561745876e-11, -0.983297266e-08,
                                   -0.265465072e-08, -0.194798427e-09, -0.513382522e-11,
                                   -0.420117906e-13, -0.469398392e-11, -0.543405219e-12,
                                   -0.121854998e-13, -0.483310746e-16, -0.429469692e-17 };

/*
    static double C1[] = {  -9e99,  0.406363373e+02,  0.291153884e+01,  0.991215929e-01,
                                    0.161603605e-02,  0.994476977e-05, -0.566497850e+01,
                                   -0.346289247e+00, -0.102486340e-01, -0.153071058e-03,
                                   -0.892381365e-06,  0.182735808e-01,  0.106282183e-02,
                                    0.311990625e-04,  0.464014079e-06,  0.269492229e-08,
                                   -0.102119482e+01, -0.649643913e-01, -0.205774955e-02,
                                   -0.323610875e-04, -0.195236396e-06,  0.531459488e-01,
                                    0.324825896e-02,  0.991819543e-04,  0.152400162e-05,
                                    0.907312536e-08, -0.132267553e-03, -0.871756401e-05,
                                   -0.262251859e-06, -0.395617938e-08, -0.232419934e-10 };

    static double C2[] = {  -9e99,  0.144323579e-02,  0.799393092e-04,  0.322526876e-05,
                                    0.596131713e-07,  0.395406097e-09, -0.839159111e-05,
                                    0.564246250e-06, -0.212045990e-07, -0.866837990e-09,
                                   -0.746255575e-11, -0.685688633e-06, -0.523054773e-07,
                                   -0.130326583e-08, -0.157964718e-10, -0.759061461e-13,
                                    0.836994324e+02, -0.609500999e+02,  0.100208335e+00,
                                   -0.688268995e+01,  0.397136599e+00, -0.250137411e-02,
                                   -0.594024621e-01,  0.457714684e-02, -0.449951913e-04,
                                   -0.273244004e+05,  0.875882129e+04, -0.227706509e+02,
                                    0.129124341e+04, -0.715722046e+02,  0.266965359e+00,
                                    0.240404391e+01, -0.269608498e+00,  0.332747493e-02 };
EQUIVALENCE (C(1),C1(1)), (C(31),C2(1))
C is both put together via thej equiv statement.
*/

    static double C[] = {  -9e99,   0.406363373e+02,  0.291153884e+01,  0.991215929e-01,
                                    0.161603605e-02,  0.994476977e-05, -0.566497850e+01,
                                   -0.346289247e+00, -0.102486340e-01, -0.153071058e-03,
                                   -0.892381365e-06,  0.182735808e-01,  0.106282183e-02,
                                    0.311990625e-04,  0.464014079e-06,  0.269492229e-08,
                                   -0.102119482e+01, -0.649643913e-01, -0.205774955e-02,
                                   -0.323610875e-04, -0.195236396e-06,  0.531459488e-01,
                                    0.324825896e-02,  0.991819543e-04,  0.152400162e-05,
                                    0.907312536e-08, -0.132267553e-03, -0.871756401e-05,
                                   -0.262251859e-06, -0.395617938e-08, -0.232419934e-10,
                                    0.144323579e-02,  0.799393092e-04,  0.322526876e-05,
                                    0.596131713e-07,  0.395406097e-09, -0.839159111e-05,
                                    0.564246250e-06, -0.212045990e-07, -0.866837990e-09,
                                   -0.746255575e-11, -0.685688633e-06, -0.523054773e-07,
                                   -0.130326583e-08, -0.157964718e-10, -0.759061461e-13,
                                    0.836994324e+02, -0.609500999e+02,  0.100208335e+00,
                                   -0.688268995e+01,  0.397136599e+00, -0.250137411e-02,
                                   -0.594024621e-01,  0.457714684e-02, -0.449951913e-04,
                                   -0.273244004e+05,  0.875882129e+04, -0.227706509e+02,
                                    0.129124341e+04, -0.715722046e+02,  0.266965359e+00,
                                    0.240404391e+01, -0.269608498e+00,  0.332747493e-02 };


    XN = XX[1]; YN = XX[2]; ZN = XX[3];
    for ( IDX=1; IDX<=5; IDX++ ) {
        X[IDX] = XN;
        Y[IDX] = YN;
        Z[IDX] = ZN;
        XN *= XX[1];
        YN *= XX[2];
        ZN *= XX[3];
    }
    g = X[1]-2.0; g2 = g*g;
    FX1 = 1.0/( 10.0 + g2 );
    XF1 = 1.0/( 15.0 - X[1] );
    g = 30.0-X[1]; g2 = g*g;
    XF2 = 1.0/g2;

    BB[1] =  Z[1]*(A[1] + A[2]*X[1] + A[3]*X[2] + A[4]*X[3] + A[5]*X[4])
               + Z[3]*(A[6] + A[7]*X[1] + A[8]*X[2] + A[9]*X[3] + A[10]*X[4])
               + Z[5]*(A[11] + A[12]*X[1] + A[13]*X[2] + A[14]*X[3] + A[15]*X[4])
               + Y[2]*Z[1]*(A[16] + A[17]*X[1] + A[18]*X[2] + A[19]*X[3] + A[20]*X[4])
               + Y[2]*Z[3]*(A[21] + A[22]*X[1] + A[23]*X[2] + A[24]*X[3] + A[25]*X[4])
               + Y[2]*Z[5]*(A[26] + A[27]*X[1] + A[28]*X[2] + A[29]*X[3] + A[30]*X[4])
               + Y[4]*Z[1]*(A[31] + A[32]*X[1] + A[33]*X[2] + A[34]*X[3] + A[35]*X[4])
               + Y[4]*Z[3]*(A[36] + A[37]*X[1] + A[38]*X[2] + A[39]*X[3] + A[40]*X[4])
               + Y[4]*Z[5]*(A[41] + A[42]*X[1] + A[43]*X[2] + A[44]*X[3] + A[45]*X[4])
               + FX1*(A[46]*Z[1] + A[47]*Z[3] + A[48]*Z[5])
               + FX1*Y[2]*(A[49]*Z[1] +  A[50]*Z[3] + A[51]*Z[5])
               + FX1*Y[4]*(A[52]*Z[1] + A[53]*Z[3] + A[54]* Z[5]);

    BB[2] = Z[1]*Y[1]*(B[1] + B[2]*X[1] + B[3]*X[2] + B[4]*X[3] + B[5]*X[4])
              + Z[3]*Y[1]*( B[6] +  B[7]*X[1] +  B[8]*X[2] +  B[9]*X[3] + B[10]*X[4])
              + Z[5]*Y[1]*(B[11] + B[12]*X[1] + B[13]*X[2] + B[14]*X[3] + B[15]*X[4])
              + Y[3]*Z[1]*(B[16] + B[17]*X[1] + B[18]*X[2] + B[19]*X[3] + B[20]*X[4])
              + Y[3]*Z[3]*(B[21] + B[22]*X[1] + B[23]*X[2] + B[24]*X[3] + B[25]*X[4])
              + Y[3]*Z[5]*(B[26] + B[27]*X[1] + B[28]*X[2] + B[29]*X[3] + B[30]*X[4])
              + Y[5]*Z[1]*(B[31] + B[32]*X[1] + B[33]*X[2] + B[34]*X[3] + B[35]*X[4])
              + Y[5]*Z[3]*(B[36] + B[37]*X[1] + B[38]*X[2] + B[39]*X[3] + B[40]*X[4])
              + Y[5]*Z[5]*(B[41] + B[42]*X[1] + B[43]*X[2] + B[44]*X[3] + B[45]*X[4]);

    BB[3] =  C[1] + C[2]*X[1] + C[3]*X[2] + C[4]*X[3] + C[5]*X[4]
                + Z[2]*( C[6] +  C[7]*X[1] +  C[8]*X[2] +  C[9]*X[3] + C[10]*X[4])
                + Z[4]*(C[11] + C[12]*X[1] + C[13]*X[2] + C[14]*X[3] + C[15]*X[4])
                + Y[2]*(C[16] + C[17]*X[1] + C[18]*X[2] + C[19]*X[3] + C[20]*X[4])
                + Y[2]*Z[2]*(C[21] + C[22]*X[1] + C[23]*X[2] + C[24]*X[3] + C[25]*X[4])
                + Y[2]*Z[4]*(C[26] + C[27]*X[1] + C[28]*X[2] + C[29]*X[3] + C[30]*X[4])
                + Y[4]*(C[31] + C[32]*X[1] + C[33]*X[2] + C[34]*X[3] + C[35]*X[4]);

    BB[3] +=  Y[4]*Z[2]*(C[36] + C[37]*X[1] + C[38]*X[2] + C[39]*X[3] + C[40]*X[4])
                + Y[4]*Z[4]*(C[41] + C[42]*X[1] + C[43]*X[2] + C[44]*X[3] + C[45]*X[4])
                + XF1*(C[46] + C[47]*Z[2] + C[48]*Z[4] + C[49]*Y[2] + C[50]*Y[2]*Z[2] + C[51]*Y[2]*Z[4] + C[52]*Y[4] + C[53]*Y[4]*Z[2] + C[54]*Y[4]*Z[4])
                + XF2*(C[55] + C[56]*Z[2] + C[57]*Z[4] + C[58]*Y[2] + C[59]*Y[2]*Z[2] + C[60]*Y[2]*Z[4] + C[61]*Y[4] + C[62]*Y[4]*Z[2] + C[63]*Y[4]*Z[4]);


    return ;

}






/*
 *  VERSION 5/13/88
 *  DEVELOPED MCDONNELL DOUGLAS
 *  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
 *
 *     INPUT -- ARGUMENT LIST
 *        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD
 *               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY
 *               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
 *               Z IS ALONG THE EARTHS NORTH DIPOLE AXIS. X IS PERPENDICULAR
 *               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
 *               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR
 *               TO X AND Z AND X Y Z FORM A RIGHT HANDED COORDINATE SYSTEM.
 *
 *     OUTPUT -- ARGUMENT LIST
 *        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE
 *               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.
 *               BB(1), BB(2) AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
 *               THE UNITS ARE NANOTESLA.
 *
 *  THIS IS THE EXPANSION FOR THE TAIL CURRENT SYSTEM. THE EXPANSION IS
 *  VALID FROM THE SUBSOLAR POIN TO -60 RE. THE EXPANSION IS BASED ON A
 *  FIT TO VALUES CALCULATED USING THE WIRE LOOP TAIL SYSTEM.
 */
void Lgm_OP88_BFTAIL( double *XX,  double *BB ) {

    double  XN, YN, ZN, R2, R, X[6], Y[6], Z[6], g, g2, R22, EXPC, TANZR, EXPR;
    int     IDX;

    static double A[] = { -9e99,  -.118386794e-12,  .260137167e+01,  .408016277e-12, -.306063863e+00,
                                   .852659791e-13,  .848404600e-14, -.568097241e-02, -.601368497e-14,
                                  -.336276159e-13, -.676779936e-15, -.110762251e-02, -.150912058e-15,
                                  -.477506548e-14, -.805245718e-02, -.130105300e-14,  .442299435e-16,
                                  -.432185140e-04, -.520612496e-01, -.918209408e-04, -.686114562e-03,
                                   .275041492e-04,  .235864029e-15, -.628394374e-04, -.236539414e-16,
                                   .379298441e-18, -.109452698e-14, -.163675727e-16, -.766199004e-04,
                                  -.110519916e-15, -.111417355e-17,  .311215382e-17, -.605957952e-06,
                                  -.609414361e+01,  .207037106e-12,  .130315144e+00, -.250115110e-13,
                                   .325228977e+00,  .169606672e-01, -.131084126e-14,  .232305257e-03,
                                   .254138418e-01, -.585580678e-03,  .344211139e-16,  .268904941e-05,
                                   .561115936e-01, -.855121118e-15,  .577135898e-03, -.389637036e-04,
                                   .531094438e-18,  .517250317e-14,  .163439821e-17,  .280008382e-15,
                                   .311491125e-17,  .165293989e-02, -.149174308e-16,  .406457779e-05,
                                  -.415855886e-06,  .127866736e-03, -.106070848e-04,  .105524883e-17,
                                   .293942950e-05, -.417367450e-06,  .134032750e-04, -.139506296e-18,
                                   0.0 };

    static double B[] = { -9e99,  -.323149328e-01,  .430535014e-02,  .115661689e-03, -.486002660e-04,
                                  -.102777234e-04, -.489864422e-05, -.356884232e-04, -.334316125e-07,
                                   .122456608e+00,  .202317315e-01, -.487990709e-03,  .338684854e-04,
                                  -.511755985e-04,  .119096933e-04,  .609353153e-03, -.243627124e-05,
                                   0.0 };

    static double C[] = {  -9e99,  .318422091e+00,  .154017442e+00,  .337581827e-01,  .436882397e-01,
                                  -.153732787e-03,  .362817457e-02,  .179382198e-03, -.394772816e-05,
                                  -.193942567e-01, -.263603775e-04, -.314364082e-04, -.103110548e-02,
                                   .386165884e-06, -.301272556e-06, -.102838611e-03, -.725608973e-04,
                                  -.893564810e-05, -.200670765e-05, -.805631807e-05, -.217861072e+02,
                                  -.219688864e+01,  .178558432e+00,  .144137907e-01, -.293171667e-04,
                                   .178727330e-01,  .846703874e-02,  .292860242e-04, -.583591628e+00,
                                   .177991433e-02,  .253212943e-02, -.629907297e-01,  .669977751e-04,
                                   .141706101e-03, -.334067698e-03,  .122648694e-03, -.259383966e-07,
                                   .252027517e-04, -.212223753e-02,  0.0 };

    XN = XX[1]; YN = XX[2]; ZN = XX[3];
    R2 = XN*XN + YN*YN + ZN*ZN;
    R = sqrt( R2 );
    for (IDX=1; IDX<=5; IDX++ ) {
        X[IDX] = XN;
        Y[IDX] = YN;
        Z[IDX] = ZN;
        XN *= XX[IDX];
        YN *= XX[IDX];
        ZN *= XX[IDX];
    }
    g = 22.0-X[1]; g2 = g*g;
    R22   = sqrt( g2 + Y[2] + Z[2] );
    EXPC  = exp( X[1]/15.0 );
    TANZR = tanh( Z[1]) * (1.0-tanh((8.0-R)/5.0) );
    g = R22-29.0; g2 = g*g;
    EXPR  = exp( -g2/60.0 );


    BB[1]=   (    A[2]*Z[1] + A[4]*X[1]*Z[1] + A[7]*Y[2]*Z[1] + A[11]*X[1]*Y[2]*Z[1]
                + A[17]*X[2]*Y[2]*Z[1] + A[18]*Z[3] + A[19]*Y[2]*Z[3] + A[20]*X[1] *Z[3]
                + A[21]*X[2]*Z[3] + A[23]*X[3]*Z[1] + A[28]*Y[4]*Z[1] + A[32]*X[4]*Z[1])*EXPC
           + (    A[33] + A[35]*X[1] + A[37]*Z[2] + A[38]*Y[2] + A[40]*Y[2]*Z[2]
                + A[41]*X[1]*Z[2] + A[42]*X[1]*Y[2] + A[44]*X[1]*Y[2]*Z[2]
                + A[45]*X[2] + A[47]*X[2]*Z[2] + A[48]*X[2]*Y[2] + A[54]*X[3]
                + A[56]*X[3] *Z[2] + A[57]*X[3]*Y[2] + A[58]*Z[4] + A[59]*Y[4]
                + A[61]*X[1]*Z[4] + A[62]*X[1]*Y[4] + A[63]*X[4])*TANZR;

    BB[2] =  (    B[1]*Y[1]*Z[1] + B[2]*X[1]*Y[1]*Z[1] + B[3]*Y[1]*Z[3] + B[4]*Y[3]*Z[1]
                + B[5]*X[1]*Y[1]*Z[3] + B[6]*X[1]*Y[3]*Z[1] + B[7]*X[2]*Y[1]*Z[1] + B[8]*X[3]*Y[1]*Z[1])*EXPC
           + (    B[9]*Y[1]*Z[1] + B[10]*X[1]*Y[1]*Z[1] + B[11]*Y[1]*Z[3] + B[12]*Y[3]*Z[1]
                + B[13]*X[1]*Y[1]*Z[3] + B[14]*X[1]*Y[3]*Z[1] + B[15]*X[2]*Y[1]*Z[1] + B[16]*X[3]*Y[1]*Z[1])*EXPR;

    BB[3] =  (    C[1] + C[2]*X[1] + C[3]*Z[2] + C[4]*Y[2] + C[5]*Y[2]*Z[2] + C[6]*X[1]*Z[2]
                + C[7]*X[1]*Y[2] + C[8]*X[1]*Y[2]*Z[2] + C[9] *X[2] + C[10]*X[2]*Z[2]
                + C[11]*X[2]*Y[2] + C[12]*X[3] + C[13]*X[3] *Z[2] + C[14]*X[3]*Y[2]
                + C[15]*Z[4] + C[16]*Y[4] + C[17]*X[1]*Z[4] + C[18]*X[1]*Y[4] + C[19]*X[4])*EXPC
           + (    C[20] + C[21]*X[1] + C[22]*Z[2] + C[23]*Y[2] + C[24]*Y[2]*Z[2] + C[25]*X[1]*Z[2] + C[26]*X[1]*Y[2]
                + C[27]*X[1]*Y[2]*Z[2] + C[28]*X[2] + C[29]*X[2]*Z[2] + C[30]*X[2]*Y[2] + C[31]*X[3] + C[32]*X[3]*Z[2]
                + C[33]*X[3]*Y[2] + C[34]*Z[4] + C[35]*Y[4] + C[36]*X[1]*Z[4] + C[37]*X[1]*Y[4] + C[38]*X[4])*EXPR;
    return;

}




/*
 *  VERSION 5/13/88
 *  DEVELOPED MCDONNELL DOUGLAS
 *  FOR INFORMATION CALL KARL PFITZER (714) 896-3231
 *
 *     INPUT -- ARGUMENT LIST
 *        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC FIELD
 *               IS TO BE DETERMINED. XX(1), XX(2), XX(3) ARE RESPECTIVELY
 *               THE X, Y, AND Z SOLAR MAGNETIC COORDINATES IN EARTH RADII.
 *               Z IS ALONG THE EARTHS NORTH DIPOLE AXIS. X IS PERPENDICULAR
 *               TO Z AND IN THE PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH
 *               LINE (X IS POSITIVE IN THE SOLAR DIRECTION). Y IS PERPENDICULAR
 *               TO X AND Z AND X,Y,Z FORM A RIGHT HANDED COORDINATE SYSTEM.
 *
 *     OUTPUT -- ARGUMENT LIST
 *        BB     A REAL ARRAY CONTAINING THE X, Y, AND Z COMPONENTS OF THE
 *               EARTH'S TOTAL MAGNETIC FIELD IN SOLAR MAGNETIC COORDINATES.
 *               BB(1), BB(2), AND BB(3) ARE THE BX, BY, AND BZ COMPONENTS.
 *               THE UNITS ARE NANOTESLA.
 *
 *  THIS SUBROUTINE CALCULATES THE FIELD FROM THE RING CURRENT SYSTEM.
 *  THE EXPANSION IS A FIT TO VALUES CALCULATED FROM THE WIRE RING
 *  CURRENT MODEL. THE EXPANSION IS VALID FROM THE SUBSOLAR POINT.
 *  TO -60 RE.
 */
void Lgm_OP88_BFRING( double *XX,  double *BB ) {

    double  XN, YN, ZN, R2, R, X[6], Y[6], Z[6], EXPC, EXPR;
    int     IDX;

    static double A[] = { -9e99,   .937029737e+00, -.734269078e+00, -.125896726e-01, -.843388063e-02,
                                   .756104711e-04,  .294507011e-02, -.719118601e-03, -.177154663e-01,
                                   .104113319e-03, -.339745485e-04,  .324439655e-03,  .492786378e-04,
                                  -.100821105e-04,  .109966887e-04,  .119616338e+00,  .403556177e+01,
                                  -.363651494e-01, -.337286459e-01, -.908902973e-04, -.980450316e-01,
                                  -.220988518e+00, -.244671475e+00, -.977974501e-03,  .311933785e-01,
                                  -.249204900e+00,  .825058070e-03,  .464195892e-02,  .223651513e-01,
                                   0.0e0 };

    static double B[] = { -9e99,  -.908641389e+00, -.249680217e-01,  .443512048e-02, -.124215709e-03,
                                   .211679921e-03, -.368134800e-04,  .547288643e-03,  .164845371e-04,
                                   .407818714e+01, -.129156231e+00, -.940633654e-01, -.220684438e+00,
                                   .878070158e-04,  .174193445e-01, -.223040987e+00,  .151981648e-01,
                                   0.0e0 };

    static double C[] = { -9e99,  -.381390073e+02, -.362173083e+01, -.410551306e+00,  .532760526e+00,
                                  -.151227645e-02,  .182345800e-01,  .358417761e-01, -.103889316e-03,
                                   .395514004e+00,  .100299786e-02,  .138275245e-03,  .288046807e-01,
                                  -.127951613e-05, -.177797800e-04,  .239511803e-02, -.284121147e-03,
                                   .939796129e-04, -.101830861e-04,  .504629929e-03,  .105982946e+02,
                                   .265464860e+01, -.157855689e+01, -.548140707e+01, -.181759612e-01,
                                   .653535097e-01,  .405331254e+00, -.726064092e-02, -.554702622e+01,
                                  -.652021402e-02,  .802389538e-01,  .167926792e+00, -.384118806e-02,
                                   .872021714e-02,  .474604567e-01,  .772720393e-01,  .144274860e-02,
                                  -.179837707e-01,  .871619151e-01,
                                   0.0e0 };

    XN = XX[1]; YN = XX[2]; ZN = XX[3];
    R2 = XN*XN + YN*YN + ZN*ZN;
    R = sqrt(R2);
    for ( IDX=1; IDX<=5; IDX++ ){
        X[IDX] = XN; Y[IDX] = YN; Z[IDX] = ZN;
        XN *= XX[1];
        YN *= XX[2];
        ZN *= XX[3];
    }
    EXPC = exp( -R/5.20 );
    if (R2 >  900.0) R2 = 900.0;
    EXPR = exp( -0.060*R2 );

    BB[1] =   (   A[ 1]*Z[ 1] + A[ 2]*X[ 1]*Z[ 1] + A[ 3]*Z[ 3] + A[ 4]*Y[ 2]*Z[ 1]
                + A[ 5]*Y[ 2]*Z[ 3] + A[ 6]*X[ 1]*Z[ 3] + A[ 7]*X[ 1]*Y[ 2]*Z[ 1]
                + A[ 8]*X[ 2]*Z[ 1] + A[ 9]*X[ 2]*Z[ 3] + A[10]*X[ 2]*Y[ 2]*Z[ 1]
                + A[11]*X[ 3]*Z[ 1] + A[12]*Z[ 5] + A[13]*Y[ 4]*Z[ 1] + A[14]*X[ 4]*Z[ 1])*EXPC
            + (   A[15]*Z[ 1] + A[16]*X[ 1]*Z[ 1] + A[17]*Z[ 3] + A[18]*Y[ 2]*Z[ 1]
                + A[19]*Y[ 2]*Z[ 3] + A[20]*X[ 1]*Z[ 3] + A[21]*X[ 1]*Y[ 2]*Z[ 1]
                + A[22]*X[ 2]*Z[ 1] + A[23]*X[ 2]*Z[ 3] + A[24]*X[ 2]*Y[ 2]*Z[ 1]
                + A[25]*X[ 3]*Z[ 1] + A[26]*Z[ 5] + A[27]*Y[ 4]*Z[ 1] + A[28]*X[ 4]*Z[ 1])*EXPR;

    BB[2] =   (   B[ 1]*Y[ 1]*Z[ 1] + B[ 2]*X[ 1]*Y[ 1]*Z[ 1] + B[ 3]*Y[ 1]*Z[ 3]
                + B[ 4]*Y[ 3]*Z[ 1] + B[ 5]*X[ 1]*Y[ 1]*Z[ 3] + B[ 6]*X[ 1]*Y[ 3]*Z[ 1]
                + B[ 7]*X[ 2]*Y[ 1]*Z[ 1] + B[ 8]*X[ 3]*Y[ 1]*Z[ 1])*EXPC
            + (   B[ 9]*Y[ 1]*Z[ 1] + B[10]*X[ 1]*Y[ 1]*Z[ 1] + B[11]*Y[ 1]*Z[ 3]
                + B[12]*Y[ 3]*Z[ 1] + B[13]*X[ 1]*Y[ 1]*Z[ 3] + B[14]*X[ 1]*Y[ 3]*Z[ 1]
                + B[15]*X[ 2]*Y[ 1]*Z[ 1] + B[16]*X[ 3]*Y[ 1]*Z[ 1])*EXPR;

    BB[3]=    (   C[ 1] + C[ 2]*X[ 1] + C[ 3]*Z[ 2] + C[ 4]*Y[ 2] + C[ 5]*Y[ 2]*Z[ 2]
                + C[ 6]*X[ 1]*Z[ 2] + C[ 7]*X[ 1]*Y[ 2] + C[ 8]*X[ 1]*Y[ 2]*Z[ 2]
                + C[ 9]*X[ 2] + C[10]*X[ 2]*Z[ 2] + C[11]*X[ 2]*Y[ 2] + C[12]*X[ 3]
                + C[13]*X[ 3]*Z[ 2] + C[14]*X[ 3]*Y[ 2] + C[15]*Z[ 4] + C[16]*Y[ 4] + C[17]*X[ 1]*Z[ 4]
                + C[18]*X[ 1]*Y[ 4] + C[19]*X[ 4])*EXPC
           + (    C[20] + C[21]*X[ 1] + C[22]*Z[ 2] + C[23]*Y[ 2] + C[24]*Y[ 2]*Z[ 2]
                + C[25]*X[ 1]*Z[ 2] + C[26]*X[ 1]*Y[ 2] + C[27]*X[ 1]*Y[ 2]*Z[ 2] + C[28]*X[ 2]
                + C[29]*X[ 2]*Z[ 2] + C[30]*X[ 2]*Y[ 2] + C[31]*X[ 3] + C[32]*X[ 3]*Z[ 2]
                + C[33]*X[ 3]*Y[ 2] + C[34]*Z[ 4] + C[35]*Y[ 4] + C[36]*X[ 1]*Z[ 4] + C[37]*X[ 1]*Y[ 4] + C[38]*X[ 4])*EXPR;

    return;

}





/*
 *  THIS FUNCTION CALCULATES THE STRENGTH OF THE RING CURRENT FROM THE
 *  STANDOFF DISTANCE AND THE DST.
 *
 *  THIS FUNCTION CAN BE USED TO CALCULATE A VALUE FOR SRING, ONE OF THE
 *  REQUIRED PARAMETERS FOR CALCUATING THE DYNAMIC MAGNETIC FIELD.
 *
 *  IT CALCULATES THE CONTRIBUTION OF THE MAGNETOPAUSE CURRENTS TO GROUND
 *  BASED SIGNATURE AND SUBTRACTS THAT COMPONENT FROM THE OBSERVED VALUE
 *  OF DST. IT ATTRIBUTES THE REMAINDER TO THE RING CURRENT.
 *
 *  INPUT PARAMETERS
 *
 *        SOFFD  THE STANDOFF DISTANCE OF THE MAGNETOPAUSE. THE QUIET STANDOFF
 *               DISTANCE IS 10.5 EARTH RADII. ACCEPTABLE VALUES RANGE BETWEEN
 *               6 AND 11. THIS VALUE IS USED TO CALCULATE THE STRENGTH OF
 *               THE MAGNETOPAUSE CURRENTS AND TO SCALE THE SIZE OF THE
 *               MAGNETOPAUSE. THIS VALUE ALSO SCALES THE SIZE OF THE TAIL
 *               CURRENT SYSTEM. THE RING SYSTEM IS NOT SCALED, SINCE ITS
 *               SOURCE IS PRIMARILY AT RADIAL DISTANCES.
 *
 *        DST    DST IS THE STANDARD PUBLISHED DST VALUE IN NANOTESLA, THE
 *               STORMTIME DISTURBED EQUATORIAL FIELD.
 *
 *  CON    SCALES THE EFFECT OF THE DST (ITS VALUE OF .03 IS STILL
 *         SOMEWHAT UNCERTAIN).
 */
double  Lgm_OP88_RINGST( double SOFFD,  double DST ) {
    double  CON, SCL, DSTMOD, RINGST;
    CON = 3.0;
    SCL = 10.5/SOFFD;
    DSTMOD = (SCL*SCL*SCL-1.0)*15.0 - DST;
    RINGST = 1.0 + DSTMOD*CON;
    return( RINGST );
}



/*
 *  THIS FUNCTION CALCULATES THE STANOFF DISTANCE FROM THE SOLAR WIND VELOCITY
 *  AND DENSITY. IT CAN BE USED TO EVALUATE THE PARAMETER SOFFD. SOFFD IS
 *  REQUIRED FOR ALL SCALING OPERATIONS.
 *
 *  INPUT PARAMETERS
 *
 *     VEL    THE SOLAR WIND VELOCITY IN KM/SEC NEAR THE SUBSOLAR POINT.
 *            TYPICAL VALUES ARE 300 TO 500.
 *
 *     DEN    THE NUMBER DENSITY OF THE SOLAR WIND IN NUMBER PER CC.
 *            TYPICAL VALUES ARE 5 TO 50.
 *
 *  OUTPUT
 *
 *     STDOFF THE DISTANCE TO THE SUBSOLAR POINT IN RE.
 *
 */
double Lgm_OP88_STDOFF( double VEL,  double DEN ) {
    double  g, g2, g3, g6, STDOFF;
    g = DEN*VEL*VEL; g2 = g*g; g3 = g*g2; g6 = g3*g3;
    STDOFF = 98.0*g6;
    return( STDOFF );
}
