import math
from ctypes import pointer, c_double

from Lgm_Wrap import Lgm_Set_Coord_Transforms, SM_TO_GSM, Lgm_Convert_Coords, TRUE, Lgm_Set_Lgm_B_OP77, Lgm_LstarInfo, SetLstarTolerances, Lgm_Trace, GSM_TO_SM, WGS84_A, RadPerDeg, NewTimeLstarInfo, Lstar
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagEphemInfo



# * Compute L*s, Is, Bms, Footprints, etc...
# * These quantities are stored in the MagEphemInfo Structure
#Lgm_Vector    Psm, P;
Psm = Lgm_Vector.Lgm_Vector(0, 0, 0)
P = Lgm_Vector.Lgm_Vector(0, 0, 0)
#Lgm_CTrans    *c = Lgm_init_ctrans(0);
c = Lgm_CTrans.Lgm_CTrans()

#Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, 1000);
MagEphemInfo = Lgm_MagEphemInfo.Lgm_MagEphemInfo(1000, 0)

#// Date and UTC
Date       = 19800625;
UTC        = 19.0;
Lgm_Set_Coord_Transforms( Date, UTC, pointer(c) );

#// Position in SM
#Psm.x = -3.0; Psm.y = 0.0; Psm.z = 0.0;
#Psm.x = -1.5; Psm.y = 0.0; Psm.z = 0.0;
#Psm.x = -1.25; Psm.y = 0.0; Psm.z = 0.0;
#Psm.x = -1.05; Psm.y = 0.0; Psm.z = 0.0;
Psm.x = -6.6
Psm.y = 0.0
Psm.z = 0.0
#Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, c );
Lgm_Convert_Coords( pointer(Psm), pointer(P), SM_TO_GSM, pointer(c) );
#// Create array of Pitch Angles to compute
#for (nAlpha=0,a=1.0; a<=90.0; a+=0.1,++nAlpha) {
#    Alpha[nAlpha] = a ;
#    printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
#}
Alpha = range(1, 90, 20)  # 1...89


#//USER INPUT STUFF
#MagEphemInfo->LstarQuality   = 3;
#MagEphemInfo->SaveShellLines = TRUE;
#MagEphemInfo->LstarInfo->VerbosityLevel = 0;
#MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;
MagEphemInfo.LstarQuality   = 3;
MagEphemInfo.SaveShellLines = TRUE;
MagEphemInfo.LstarInfo.contents.VerbosityLevel = 0;
MagEphemInfo.LstarInfo.contents.mInfo.contents.VerbosityLevel = 0;


Kp = 1;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
#MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
Lgm_Set_Lgm_B_OP77( MagEphemInfo.LstarInfo.contents.mInfo )
#MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
MagEphemInfo.LstarInfo.contents.mInfo.contents.Kp = Kp
#if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;



#Lgm_LstarInfo 	*LstarInfo, *LstarInfo2, *LstarInfo3;
LstarInfo = Lgm_LstarInfo()
#Lgm_Vector  v, v1, v2, v3, vv1, Bvec;

#LstarInfo = MagEphemInfo->LstarInfo;
LstarInfo = MagEphemInfo.LstarInfo;


#// Save Date, UTC to MagEphemInfo structure
#MagEphemInfo->Date   = Date;
#MagEphemInfo->UTC    = UTC;
MagEphemInfo.Date   = Date
MagEphemInfo.UTC    = UTC

#// Save nAlpha, and Alpha array to MagEphemInfo structure
#MagEphemInfo->nAlpha = nAlpha;
#for (i=0; i<MagEphemInfo->nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];
MagEphemInfo.nAlpha = len(Alpha)

#for (i=0; i<MagEphemInfo->nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];
MagEphemInfo.Alpha = (c_double*len(Alpha))(*Alpha)


#// Set Tolerances
#SetLstarTolerances( MagEphemInfo->LstarQuality, LstarInfo );
SetLstarTolerances(MagEphemInfo.LstarQuality, MagEphemInfo.LstarInfo )


#sprintf( Filename, "DipoleTest_6.6/results_%.0e.dat", LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
#fpout = fopen(Filename, "w");
print("DipoleTest_6.6/results_%.0e.dat" %
      (MagEphemInfo.LstarInfo.contents.mInfo.contents.Lgm_FindShellLine_I_Tol) )


#// set coord transformation
#Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );
Lgm_Set_Coord_Transforms(Date, UTC, MagEphemInfo.LstarInfo.contents.mInfo.contents.c)

# *  Blocal at sat location.
#MagEphemInfo->P = *u;
MagEphemInfo.P = P

Bvec = Lgm_Vector.Lgm_Vector(0,0,0)
#LstarInfo->mInfo->Bfield( u, &Bvec, LstarInfo->mInfo );
MagEphemInfo.LstarInfo.contents.mInfo.contents.Bfield(pointer(P), pointer(Bvec),
                                                      MagEphemInfo.LstarInfo.contents.mInfo)
#Blocal = Lgm_Magnitude( &Bvec );
Blocal = Bvec.magnitude()

#MagEphemInfo->B = Blocal;
MagEphemInfo.B = Blocal

u=P
v1 = Lgm_Vector.Lgm_Vector(0,0,0)
v2 = Lgm_Vector.Lgm_Vector(0,0,0)
v3 = Lgm_Vector.Lgm_Vector(0,0,0)
vv1 = Lgm_Vector.Lgm_Vector(0,0,0)
TRACE_TOL = 1e-7

# *  Compute Field-related quantities for each Pitch Angle.
#if ( Lgm_Trace( u, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo->mInfo ) == 1 ) {

if Lgm_Trace(pointer(u), pointer(v1), pointer(v2), pointer(v3),
          120, 0.01,
          TRACE_TOL, MagEphemInfo.LstarInfo.contents.mInfo) == 1:
    #    MagEphemInfo->Pmin = v3;
    MagEphemInfo.Pmin = v3
#    MagEphemInfo->Bmin     = LstarInfo->mInfo->Bmin;
    MagEphemInfo.Bmin = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bmin
#         *  Get a simple measure of how big L is
#        Lgm_Convert_Coords( &v1, &vv1, GSM_TO_SM, LstarInfo->mInfo->c );
    Lgm_Convert_Coords( pointer(v1), pointer(vv1), GSM_TO_SM,
                       MagEphemInfo.LstarInfo.contents.mInfo.contents.c );

#        Lam = asin( vv1.z/Lgm_Magnitude( &vv1 ) );
    Lam = math.asin(vv1.z/vv1.magnitude())
#        CosLam = cos( Lam );
    CosLam = math.cos(Lam)
#        LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam );
    LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam )

#        for ( i=0; i<MagEphemInfo->nAlpha; i++ ){  // LOOP OVER PITCH ANGLES
    for i, pa in enumerate(Alpha):
#            // colorize the diagnostic messages.
#            sprintf( LstarInfo3->PreStr, "\033[38;5;%dm", Colors[i%9]); sprintf( LstarInfo3->PostStr, "\033[0m");
        print(MagEphemInfo.LstarInfo.contents.PreStr)
        print(MagEphemInfo.LstarInfo.contents.PostStr)
#            PreStr = LstarInfo3->PreStr; PostStr = LstarInfo3->PostStr;
        PreStr = MagEphemInfo.LstarInfo.contents.PreStr
        PostStr = MagEphemInfo.LstarInfo.contents.PostStr

#             *  Set Pitch Angle, sin, sin^2, and Bmirror
#            sa = sin( MagEphemInfo->Alpha[i]*RadPerDeg ); sa2 = sa*sa;
        sa = math.sin( pa*RadPerDeg )
        sa2 = sa*sa;

#            printf("%sComputing L* for Pitch Angle: Alpha[%d] = %g Date: %ld   UTC: %g   Lsimple = %g%s\n", PreStr, i, MagEphemInfo->Alpha[i], Date, UTC, LSimple, PostStr );
        print("%sComputing L* for Pitch Angle: Alpha[%d] = %g Date: %d   UTC: %g   Lsimple = %g%s\n" %
              ( PreStr, i, MagEphemInfo.Alpha[i], Date, UTC, LSimple, PostStr ))

#            LstarInfo3->mInfo->Bm = Blocal/sa2;
        MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm = Blocal/sa2
#            NewTimeLstarInfo( Date, UTC, MagEphemInfo->Alpha[i], LstarInfo3->mInfo->Bfield, LstarInfo3 );
        Lgm_Set_Coord_Transforms( Date, UTC, MagEphemInfo.LstarInfo.contents.mInfo.contents.c )
        MagEphemInfo.LstarInfo.contents.PitchAngle = pa

#            MagEphemInfo->Bm[i] = LstarInfo3->mInfo->Bm;
        MagEphemInfo.Bm[i] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm

#             *  Compute L*

#            if ( LSimple < 10.0 ){
        if LSimple < 10.0:
#// USER SHOULD DECIDE THRESHOLD HERE
#
#                LstarInfo2->mInfo->Bm = LstarInfo3->mInfo->Bm;
#                if (LstarInfo3->VerbosityLevel >= 2 ) {
            if MagEphemInfo.LstarInfo.contents.VerbosityLevel >= 2 :
#                    printf("\n\n\t%sComputing L* for: UTC = %g PA = %d  (%g)%s\n", PreStr, UTC, i, MagEphemInfo->Alpha[i], PostStr );
#                    printf("    \t%s                  I   = %g PA = %d  (%g)%s\n", PreStr, MagEphemInfo->I[i], i, MagEphemInfo->Alpha[i], PostStr );

                print("\n\n\t%sComputing L* for: UTC = %g PA = %d  (%g)%s\n" %
                      ( PreStr, UTC, i, MagEphemInfo.Alpha[i], PostStr ))
                print("    \t%s                  I   = %g PA = %d  (%g)%s\n" %
                    (PreStr, MagEphemInfo.I[i], i, MagEphemInfo.Alpha[i], PostStr ))
#                }
#                LS_Flag = Lstar( &v3, LstarInfo2);
            LS_Flag = Lstar( pointer(v3), MagEphemInfo.LstarInfo)

#                if (LstarInfo3->VerbosityLevel >= 2 ) {
#                    printf("\t%sUTC, L*          = %g %g%s\n", PreStr, UTC, LstarInfo2->LS, PostStr );
#                    printf("\t%sUTC, L*_McIlwain = %g %g%s\n", PreStr, UTC, LstarInfo2->LS_McIlwain_M, PostStr );
#                    printf("\t%sUTC, LSimple     = %g %g%s\n\n\n", PreStr, UTC, LSimple, PostStr );
#                }
            if MagEphemInfo.LstarInfo.contents.VerbosityLevel >= 2:
                print("\t%sUTC, L*          = %g %g%s\n" %
                       (PreStr, UTC, MagEphemInfo.LstarInfo.contents.LS, PostStr ))
                print("\t%sUTC, L*_McIlwain = %g %g%s\n" %
                       (PreStr, UTC, MagEphemInfo.LstarInfo.contents.LS_McIlwain_M, PostStr ))
                print("\t%sUTC, LSimple     = %g %g%s\n\n\n" %
                       (PreStr, UTC, LSimple, PostStr ))

#                MagEphemInfo->Lstar[i] = LstarInfo2->LS;
            MagEphemInfo.Lstar[i] = MagEphemInfo.LstarInfo.contents.LS

#
#
#                /*
#                 * Save results to the MagEphemInfo structure.
#                 */
#                MagEphemInfo->nShellPoints[i] = LstarInfo2->nPnts;
            MagEphemInfo.nShellPoints[i] = MagEphemInfo.LstarInfo.contents.nPnts
#                for (nn=0; nn<LstarInfo2->nPnts; nn++ ){
            for nn in range(MagEphemInfo.LstarInfo.contents.nPnts):
#                    MagEphemInfo->ShellI[i][nn] = LstarInfo2->I[nn];
                MagEphemInfo.ShellI[i][nn] = MagEphemInfo.LstarInfo.contents.I[nn]
#                    MagEphemInfo->ShellEllipsoidFootprint_Pn[i][nn] = LstarInfo2->Ellipsoid_Footprint_Pn[nn];
                MagEphemInfo.ShellEllipsoidFootprint_Pn[i][nn] = MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Pn[nn]
#                    MagEphemInfo->ShellEllipsoidFootprint_Ps[i][nn] = LstarInfo2->Ellipsoid_Footprint_Ps[nn];
                MagEphemInfo.ShellEllipsoidFootprint_Ps[i][nn] = MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Ps[nn]
#                    MagEphemInfo->ShellMirror_Pn[i][nn]    = LstarInfo2->Mirror_Pn[nn];
                MagEphemInfo.ShellMirror_Pn[i][nn] = MagEphemInfo.LstarInfo.contents.Mirror_Pn[nn]
#                    MagEphemInfo->ShellMirror_Ps[i][nn]    = LstarInfo2->Mirror_Ps[nn];
                MagEphemInfo.ShellMirror_Ps[i][nn] = MagEphemInfo.LstarInfo.contents.Mirror_Ps[nn]
#                    MagEphemInfo->ShellMirror_Ss[i][nn]    = LstarInfo2->mInfo->Sm_South;
                MagEphemInfo.ShellMirror_Ss[i][nn] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Sm_South
#                    MagEphemInfo->ShellMirror_Sn[i][nn]    = LstarInfo2->mInfo->Sm_North;
                MagEphemInfo.ShellMirror_Sn[i][nn] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Sm_North

#
#                    /*
#                     *  Save all of the drift shell FLs in MagEphemInfo structure
#                     */
#                    MagEphemInfo->nFieldPnts[i][nn] = LstarInfo2->nFieldPnts[nn];
                MagEphemInfo.nFieldPnts[i][nn] = MagEphemInfo.LstarInfo.contents.nFieldPnts[nn]
#                    for (tk=0; tk<LstarInfo2->nFieldPnts[nn]; tk++){    // loop over points in a FL
                for tk in range(MagEphemInfo.LstarInfo.contents.nFieldPnts[nn]):
#                        MagEphemInfo->s_gsm[i][nn][tk] = LstarInfo2->s_gsm[nn][tk];
                    MagEphemInfo.s_gsm[i][nn][tk] = MagEphemInfo.LstarInfo.contents.s_gsm[nn][tk]
#                        MagEphemInfo->Bmag[i][nn][tk]  = LstarInfo2->Bmag[nn][tk];
                    MagEphemInfo.Bmag[i][nn][tk] = MagEphemInfo.LstarInfo.contents.Bmag[nn][tk]
#                        MagEphemInfo->x_gsm[i][nn][tk] = LstarInfo2->x_gsm[nn][tk];
                    MagEphemInfo.x_gsm[i][nn][tk] = MagEphemInfo.LstarInfo.contents.x_gsm[nn][tk]
#                        MagEphemInfo->y_gsm[i][nn][tk] = LstarInfo2->y_gsm[nn][tk];
                    MagEphemInfo.y_gsm[i][nn][tk] = MagEphemInfo.LstarInfo.contents.y_gsm[nn][tk]
#                        MagEphemInfo->z_gsm[i][nn][tk] = LstarInfo2->z_gsm[nn][tk];
                    MagEphemInfo.z_gsm[i][nn][tk] = MagEphemInfo.LstarInfo.contents.z_gsm[nn][tk]
#                    }
#                }
        else:
            print(" Lsimple >= 10.0  ( Not doing L* calculation )" )
#            } else {
#                printf(" Lsimple >= 10.0  ( Not doing L* calculation )\n" );
#            }
#fprintf(fpout, "%.15lf %.15lf\n", MagEphemInfo->Alpha[i], (6.6-LstarInfo2->LS)/LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
    print("%.15lf %.15lf\n" %
          ( MagEphemInfo.Alpha[i], (6.6-MagEphemInfo.LstarInfo.contents.LS/MagEphemInfo.LstarInfo.contents.mInfo.contents.Lgm_FindShellLine_I_Tol)))

print(MagEphemInfo.P.x, MagEphemInfo.P.y, MagEphemInfo.P.z)
for i in range(len(Alpha)):
    print(Alpha[i], MagEphemInfo.Lstar[i])
