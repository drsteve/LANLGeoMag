#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <Lgm_QinDenton.h>
#include <fcntl.h>


int main(){
    long int            Date;
    double              L, I, Bm, M, a, UTC, JD;
    double              r, lat, lon;
    Lgm_Vector          u, v;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_QinDentonOne    p;
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );

    mInfo->Bfield = Lgm_B_TS04;
    mInfo->P      = 4.1011111111111118;
    mInfo->Dst    = 7.7777777777777777;
    mInfo->By     = 3.7244444444444444;
    mInfo->Bz     = -0.12666666666666665;
    mInfo->W[0]   = 0.12244444444444445;
    mInfo->W[1]   = 0.2514;
    mInfo->W[2]   = 0.089266666666666661;
    mInfo->W[3]   = 0.047866666666666668;
    mInfo->W[4]   = 0.22586666666666666;
    mInfo->W[5]   = 1.0461333333333334;

    //vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
    Date = 19960106;
    UTC  = 1.2444444444444445;
Date = 20130101;
UTC = 0.0;
Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );    // Compute JD
// Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
Lgm_get_QinDenton_at_JD( JD, &p, 0 );
// Set params in mInfo structure.
Lgm_set_QinDenton( &p, mInfo );


mInfo->Bfield = Lgm_B_T89;
mInfo->Kp = 3;




    Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    r = 4.83415065;
    lat = 36.44696388*RadPerDeg;
    lon = -40.26632902*RadPerDeg;
    u.x = r*cos(lat)*cos(lon);
    u.y = r*cos(lat)*sin(lon);
    u.z = r*sin(lat);
//u.x = 0.3503119221132272;
//u.y = 0.185820103265288;
//u.z = -1.0028377093930358;
    printf("u_gsm = %.15lf %.15lf %.15lf\n", u.x, u.y, u.z);
//    Lgm_Convert_Coords( &u, &v, GSM_TO_WGS84, c );
//    printf("v_wgs84 = %.15lf %.15lf %.15lf\n", v.x, v.y, v.z);
    a = 90.0;

    L = Lgm_McIlwain_L( Date, UTC, &u, a, 1, &I, &Bm, &M, mInfo );
    printf("Pitch Angle: %g    McIlwain L  = %.15g   ( I, Bm, M = %.15g %g %g )\n", a, L, I, Bm, M);


    Lgm_FreeMagInfo( mInfo );
    Lgm_free_ctrans( c );



    return(0);
}


