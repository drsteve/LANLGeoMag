/* Clang port of TA16 model with LANLGeoMag hooks
 * Originally ported from F77 to F2008 by S.K. Morley, Nov. 2021
 * This port by S.K. Morley, July 2022
 * A data-based RBF model driven by interplanetary and ground-based data
 * Ported from version dated 10/19/2016, including the correction of Nov. 2021
 * which was found during the porting process.
 * This model is based on fitting to a data set with 732,746 records
 * References: (1) Tsyganenko, N. A., and V. A. Andreeva (2016), "An empirical RBF model of the magnetosphere
 *                 parameterized by interplanetary and ground-based drivers", v.121, doi:10.1002/2016JA023217,
 *                 accepted by JGRA, 10/17/2016.
 *             (2) Andreeva, V. A., and N. A. Tsyganenko (2016), "Reconstructing the magnetosphere from data 
 *                 using radial basis functions, JGRA Space Physics, v.121, 2249-2263, doi:10.1002/2015JA022242.
 * ORIGINAL CODE BY: N. A. Tsyganenko AND V. A. Andreeva
 *------------------------------------------------------------------------------------------------------
 * NOTE: THE MODEL IS VALID ONLY UP TO Xgsw=-15 Re and should NOT be used tailward of that distance
 *------------------------------------------------------------------------------------------------------
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_QinDenton.h"
#include "Lgm/Lgm_TA2016.h"
#ifndef LGM_INDEX_DATA_DIR
#warning "hard-coding LGM_INDEX_DATA_DIR because it was not in config.h"
#define LGM_INDEX_DATA_DIR    /usr/local/share/LanlGeoMag/Data
#endif

static double Lgm_TA16_As[22] = {
    12.544, -0.194, 0.305, 0.0573, 2.178,
    0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0, 0.0, 
    -0.210, 0.0, 0.0, -0.636, -2.600, 0.832, -5.328, 
    1.103, -0.907, 1.450};


void Lgm_Init_TA16(LgmTA16_Info *ta){
    ta->SetupDone = FALSE;
    ta->readTAparams = TRUE;  //TRUE to read from Tsyganenko's files
                              // FALSE to estimate from QD data
    ta->lenA = 23329;  //23328 parameter values, 1-based indexing
    ta->lenL = 1297;  // max val of L in fortran code is 1296
    LGM_ARRAY_1D( ta->A,   ta->lenA, double );
    LGM_ARRAY_1D( ta->XX,   ta->lenL, double );
    LGM_ARRAY_1D( ta->YY,   ta->lenL, double );
    LGM_ARRAY_1D( ta->ZZ,   ta->lenL, double );
    LGM_ARRAY_1D( ta->ST,   ta->lenL, double );
    LGM_ARRAY_1D( ta->RHO,   ta->lenL, double );
    LGM_ARRAY_1D( ta->ZSP,   ta->lenL, double );
    LGM_ARRAY_1D( ta->ZCP,   ta->lenL, double );
    LGM_ARRAY_1D( ta->RHBR,   ta->lenL, double );
    ta->ArraysAlloced = TRUE;
    Lgm_GetData_TA16(ta);
    TA2016_SetGrid(ta);
    return;
}

void Lgm_DeAllocate_TA16( LgmTA16_Info *t ){

    if ( t->ArraysAlloced == TRUE ) {
        LGM_ARRAY_1D_FREE( t->A );
        LGM_ARRAY_1D_FREE( t->XX );
        LGM_ARRAY_1D_FREE( t->YY );
        LGM_ARRAY_1D_FREE( t->ZZ );
        LGM_ARRAY_1D_FREE( t->ST );
        LGM_ARRAY_1D_FREE( t->RHO );
        LGM_ARRAY_1D_FREE( t->ZSP );
        LGM_ARRAY_1D_FREE( t->ZCP );
        LGM_ARRAY_1D_FREE( t->RHBR );
        t->ArraysAlloced = FALSE;
    }
}

int Lgm_Copy_TA16_Info(LgmTA16_Info *targ, LgmTA16_Info *src) {
    int i;
    if ( src == NULL) {
        printf("Lgm_Copy_TA16_Info: Error, source structure is NULL\n");
        return(-1);
    }
    memcpy(targ, src, sizeof(LgmTA16_Info));
    LGM_ARRAY_1D( targ->A,   targ->lenA, double );
    LGM_ARRAY_1D( targ->XX,  targ->lenL, double );
    LGM_ARRAY_1D( targ->YY,  targ->lenL, double );
    LGM_ARRAY_1D( targ->ZZ,  targ->lenL, double );
    LGM_ARRAY_1D( targ->ST,  targ->lenL, double );
    LGM_ARRAY_1D( targ->RHO, targ->lenL, double );
    LGM_ARRAY_1D( targ->ZSP, targ->lenL, double );
    LGM_ARRAY_1D( targ->ZCP, targ->lenL, double );
    LGM_ARRAY_1D( targ->RHBR,targ->lenL, double );
    targ->ArraysAlloced = TRUE;
    for (i=0; i<targ->lenA+1; i++) {
      targ->A[i] = src->A[i];
    }
    for (i=0; i<targ->lenL+1; i++) {
      targ->XX[i] = src->XX[i];
      targ->YY[i] = src->YY[i];
      targ->ZZ[i] = src->ZZ[i];
      targ->ST[i] = src->ST[i];
      targ->RHO[i] = src->RHO[i];
      targ->ZSP[i] = src->ZSP[i];
      targ->ZCP[i] = src->ZCP[i];
      targ->RHBR[i] = src->RHBR[i];
    }
    return(1);
}


int Lgm_GetData_TA16(LgmTA16_Info *ta) {
    FILE *fp; 
    char line[64], Filename[1024];
    int idx=1, maxidx=23329;

    sprintf(Filename, "%s/TA_DATA/TA16_RBF.par", LGM_INDEX_DATA_DIR);
    if ((fp = fopen(Filename, "r")) != NULL) {
        for (idx; idx<maxidx; ++idx) {
            fscanf(fp, "%lf", &ta->A[idx]);
        }
    }
    fclose(fp);
    return(1);
}


int Lgm_SetCoeffs_TA16(long int Date, double UTC, LgmTA16_Info *ta) {
    double symh, bt, bz2, by2;
    double *symhc, *newell_norm, *clock;
    int idx, ncols;
    FILE *fp;
    char *Path, TA16_path[2048], datafile[2048], tmpstr[1300];
    int inyear, inmonth, inday, indoy;
    int year, month, day, doy, hour, minute;
    double bx_av, by_av, bz_av, vx, vy, vz, nden, temp;
    int IMFflag, SWflag;
    double tilt, pdyn, newe_avg, boyn_avg, symhc_avg;
    double lineutc;
    double fivemin=5./60.0-2.78e-10;  //five minutes minus a microsecond
    if ( !(ta->SetupDone) ){
        Lgm_Init_TA16( ta );
    }
    Lgm_Doy( Date, &inyear, &inmonth, &inday, &indoy );

    // Set default values
    ta->Pdyn = 3.;  // Instantaneous Pdyn
    symh = -29.354176;  // when corrected gives -46nT (6dp)
    // Sym-Hc is pressure corrected. Input is 30-minute
    // avg. of Sym-Hc centered on observation time
    ta->SymHc_avg = 0.8 * symh - 13*sqrt(ta->Pdyn);
    ta->Xind_avg = 1.;  //Avg. of coupling param based on Newell, over previous 30 minutes
    ta->By_avg = 3.;  // Avg. over previous 30 minutes of By
    // TODO: read model parameters from file.

    // If flag set, calculate from Qin-Denton...
    if ( !(ta->readTAparams) ){
      Lgm_QinDenton *qin=Lgm_init_QinDenton(0);
      Lgm_read_QinDenton(Date, qin);
      // Init array for pressure corrected sym-h
      symhc = (double *)calloc( qin->nPnts, sizeof(double) );
      newell_norm = (double *)calloc( qin->nPnts, sizeof(double) );
      for (idx=0; idx<qin->nPnts; idx++) {
          bz2 = qin->BzIMF[idx]*qin->BzIMF[idx];
          by2 = qin->ByIMF[idx]*qin->ByIMF[idx];
          bt = sqrt(bz2 + by2);  //QD files don't have Bx, so just use transverse
          symhc[idx] = 0.8 * qin->Dst[idx] - 13*sqrt(qin->Pdyn[idx]);
          newell_norm[idx] = 1.0e-4 * pow(qin->V_SW[idx], 4.0/3.0) *
                                      pow(bt, 2.0/3.0) *
                                      pow(sin(clock[idx]/2.0), 8.0/3.0);
          //TODO: take window averages for given time
      }
      // Clean up
      Lgm_destroy_QinDenton(qin);
      free(symhc);
      free(newell_norm);
    } else {
      // check that data location is set and load TA annual file
      // use nearest point...
      Path = getenv("TA16_DATA_PATH");
      if (Path ==NULL) {  // setting path
          strcpy( TA16_path, LGM_INDEX_DATA_DIR );
          strcat( TA16_path, "/TA_DATA" );
      } else {
        if ( access( TA16_path, F_OK ) < 0 ) {
            printf("Lgm_SetCoeffs_TA16: Warning, TA16 Data directory not found at %s. Use TA16_DATA_PATH environment variable to set path.\n", TA16_path );
            exit(-1);
        } else{
          strcpy( TA16_path, Path );
        }
      }  // done setting path

      // Now load data filea
      // IYEAR       i4       4-digit year
      // IDAY        i4       Day of year (IDAY=1 is January 1)
      // IHOUR       i3       UT hour (0 - 23)
      // MIN         i3       UT minute (0 - 59)
      // <Bx IMF>    f8.2     nT, in GSW coordinates, average over 30-min trailing interval
      // <By IMF>    f8.2     nT, in GSW coordinates, average over 30-min trailing interval
      // <Bz IMF>    f8.2     nT, in GSW coordinates, average over 30-min trailing interval
      // Vx          f8.1     solar wind Vx in GSE coordinates, km/s
      // Vy          f8.1     solar wind Vy in GSE coordinates, km/s
      // Vz          f8.1     solar wind Vz in GSE coordinates, km/s
      // Np          f7.2     solar wind proton density, cm^-3
      // T           f9.0     solar wind proton temperature, degs K
      // Sym-H       f7.1     Sym-H index, nT
      // IMF flag    i5       equals 1 for actually measured or 2 for linearly interpolated IMF
      // SW flag     i5       equals 1 for actually measured or 2 for linearly interpolated solar wind parameters
      // Tilt        f8.4     geodipole tilt angle (radians!) in GSW coordinate system
      // Pdyn        f7.2     solar wind flow ram pressure, nPa
      // <N-index>   f8.4     average over 30-min trailing interval, see Eq.(1) in TA16_Model_description.pdf 
      // <B-index>   f8.4     average over 30-min trailing interval, see Eq.(2) in TA16_Model_description.pdf
      // <SymHc>     f7.1     sliding average over 30-min interval, centered on the current time moment
      sprintf( datafile, "%s/%ld_OMNI_5m_with_RBF_TA16_drivers.dat", TA16_path, Date/10000);
      if ( (fp = fopen( datafile, "r" )) != NULL ) {
        // to start with, just loop over... should we actually
        // be interpolating linearly between values?
        while ( fgets( &tmpstr, 1300, fp ) != NULL ) {
            ncols = sscanf( tmpstr, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf",
                            &year, &doy, &hour, &minute, &bx_av, &by_av, &bz_av,
                            &vx, &vy, &vz, &nden, &temp, &symh, &IMFflag, &SWflag,
                            &tilt, &pdyn, &newe_avg, &boyn_avg, &symhc_avg);
            lineutc = (double)hour + (double)minute/60.0;
            if ( (abs(doy - indoy) == 0) && (fabs(lineutc - UTC) < fivemin) ) {
                //nearest time, so fill ta structure
                printf("Using %s", tmpstr);
                ta->Pdyn = pdyn;
                ta->Xind_avg = newe_avg;
                ta->By_avg = by_av;
                ta->SymHc_avg = symhc_avg;
                break;
            }
        }
      } else {//end IF for file opened successfully
        printf("Lgm_SetCoeffs_TA16(): Line %d in file %s. Could not open file %s\n", __LINE__, __FILE__, datafile );
        exit(-1);
      }
    }  // done loading TA files
    // TODO: make sure I have fallbacks to default values

}


int TA2016_SetGrid(LgmTA16_Info *Info) {
    // initialize the TA2016 coefficient set, grid, etc.
    int Klat, Nlat, Nlon, I, J, K, L=0;
    double D=4.0;
    double Alpha, XXXX, YYYY, ZZZZ, RHighGrid, RLowGrid, X4sq, Y4sq, Z4sq;
    double Xlon, XlonD, Xcolat, XcolatD, sinXco, sinXlo, cosXlo;
    double R0, R1, RH, RM, Rho1, R, RLAST, STS, STT, STN;
    double CN, CS, CTN, CTS, CTT, DN, DS, EN, ES, CP, SP;
    double T, PSIN, PSIS, dLonDeg;
    double B0, F1, PD_TR, PM=0.0;
  
    Klat  = 8;          // Klat IS THE NUMBER OF LATITUDE CIRCLES IN THE NORTHERN HEMISPHERE (EXCLUDING THE POLE)
    RLowGrid  = 3.3;  //  THE INNERMOST SPHERE RADIUS
    RHighGrid =16.0;   //  UPPER LIMIT ON THE RADIUS OF OUTERMOST RBF SPHERE
    // Upper limit on the data tailward limit = 20.0
  
    RH=8.0;        // HINGING PARAMETERS (SEE TAG-2015)
    Alpha=3.0;     // Alpha PARAMETER (SEE TAG-2015)
    Nlat = Klat+1;   // NUMBER OF LATITUDE CIRCLES IN THE N. HEMISPHERE INCLUDING THE POLE
  
    R = RLowGrid;  //  THE INNERMOST SPHERE RADIUS
  
    PD_TR = 0.5;    //  this is just to filter out RBF centers outside the Lin et al model magnetopause
    CTN = cos(Lgm_TA16_As[19]);
    CTS = CTN;
    STN = sin(Lgm_TA16_As[19]);
    STS = STN;
    CN = 0.0;
    CS = 0.0;
  
    for (J=1; J<=100; J++) {   // J COUNTS THE NUMBER OF SPHERES WITH RBF CENTERS
                 // (WHICH IS ACTUALLY MUCH LESS THAN 100)
      for (I=1; I<=Nlat; I++) {     // I COUNTS THE LATITUDE CIRCLES (FROM NORTH POLE DOWN TO EQUATOR)
        XcolatD = (90.0/(Nlat-0.5))*((double)(I)-1.0);  //  COLATITUDE OF THE Ith CIRCLE (in degs)
        Nlon=4*(I-1);                       //  NUMBER OF RBF CENTERS ON THE ItH CIRCLE
        if (I != 1) {
          dLonDeg = 360.0/Nlon;            //  LONGITUDE INTERVAL BETWEEN CENTERS (degs)
        } else {
          Nlon = 1;                         //  NUMBER OF RBF CENTERS ON THE NORTH POLE
          dLonDeg = 0.0;
        }
        for (K=1; K<=Nlon; K++) {
          XlonD = (K-1)*dLonDeg;
          // number below is radian conversion... do with named variable
          Xcolat = XcolatD*0.01745329252;  // COLATITUDE and
          Xlon = XlonD*0.01745329252;      // LONGITUDE OF THE Lth RBF CENTER

          sinXco = sin(Xcolat);
          sinXlo = sin(Xlon);
          cosXlo = cos(Xlon);
          XXXX = R*sinXco*cosXlo;
          YYYY = R*sinXco*sinXlo;
          ZZZZ = R*cos(Xcolat);
          // HERE WE CALCULATE THE DISTANCE FROM THE RBF NODE TO THE LIN ET AL. [2010] MODEL MAGNETOPAUSE.
          // THE RBF NODES LYING OUTSIDE THE BOUNDARY ARE NOT INCLUDED IN THE GRID
          EN = Lgm_TA16_As[21];
          ES = Lgm_TA16_As[21];
  
          X4sq = XXXX*XXXX;
          Y4sq = YYYY*YYYY;
          Z4sq = ZZZZ*ZZZZ;
          R1 = sqrt(X4sq+Y4sq+Z4sq);
          Rho1 = sqrt(Y4sq+Z4sq);

          CTT=XXXX/R1;
          STT = sqrt(Y4sq+Z4sq)/R1;
          T = atan2(STT, CTT);
          if (Rho1 > 1.0e-5) {  // WE ARE NOT ON THE X-AXIS - NO SINGULARITIES TO WORRY ABOUT
            SP = ZZZZ/Rho1;
            CP = YYYY/Rho1;
          } else {              // WE ARE ON THE X-AXIS
            if (XXXX > 0.0) {   // IF ON THE DAYSIDE, NO PROBLEM EITHER
              SP = 0.0;
              CP = 1.0;
            } else {            // WE ARE ON THE TAIL AXIS; TO AVOID SINGULARITY:
              RM = 1000.0;      // ASSIGN RM=1000 (A CONVENTIONAL SUBSTITUTE VALUE)
              // Raise an error here? Or set bad values and a warning?
              return(-1);       //  AND EXIT
            }
          }
  
          PSIN = acos(CTT*CTN+STT*STN*SP);
          PSIS = acos(CTT*CTS-STT*STS*SP);
  
          DN = Lgm_TA16_As[16];
          DS = Lgm_TA16_As[16];
  
          B0 = Lgm_TA16_As[6];
          F1 = pow(sqrt(0.5*(1.0+CTT))+Lgm_TA16_As[5]*2.0*STT*CTT*(1.0-exp(-T)), B0);
          R0 = Lgm_TA16_As[0]*pow(PD_TR+PM, Lgm_TA16_As[1]);
          RM = R0*F1+CN*exp(DN*pow(PSIN, EN))+CS*exp(DS*pow(PSIS, ES));    // POSITION OF THE MODEL MAGNETOPAUSE
          if (R <= RM) {
            L = L+1;              //  COUNTER OF THE RBF CENTERS
            Info->XX[L] = XXXX;
            Info->YY[L] = YYYY;
            Info->ZZ[L] = ZZZZ;
            Info->ST[L] = sinXco;
            Info->RHO[L] = R*sinXco;
            Info->ZSP[L] = ZZZZ*sinXlo;
            Info->ZCP[L] = ZZZZ*cosXlo;
            Info->RHBR[L] = RH/R*(1.0 - pow(1.0 + pow(R/RH, Alpha), 1.0/Alpha));
          }
        }  // End loop over K
      }  // End loop over I
  
      RLAST = R;
      R = R*(Nlat-0.5+M_PI/4.0)/(Nlat-0.5-M_PI/4.0);  // Increment R by a fixed factor
      if (R > RHighGrid) break;                       // RBF centers only inside R=RHIGH
    }  // End loop over J
  Info->SetupDone = TRUE ;
}


int TA2016(Lgm_Vector *posGSM, double *PARMOD, Lgm_CTrans *ctrans, Lgm_Vector *BvecGSM, LgmTA16_Info *Info) {
    int I;
    Lgm_Vector posSM, BvecSM;
    Lgm_Vector DP, DM, Parr, Marr, PCP, PCM, DCM, DCP, DCMsq, DCPsq;
    Lgm_Vector TCM, TCP;
    Lgm_Vector CTarr, CParr, STarr, SParr;
    double diptilt, PS, cPS, sPS, tPS, CM, CP;
    double Pdyn, SymV, Xind, ByIMF, FPD;
    double ACP, ACT, AP, ASP,AST, AT, sDT, cDTM1;
    double DCMXY, DCMXZ, DCMYZ, DCPXY, DCPXZ, DCPYZ, DELTA_ZR, DTHETA;
    double D=4.0;
    int x=1, y=2, z=3;

    diptilt = ctrans->psi;
    cPS = cos(diptilt);  // cos(dipole tilt)
    sPS = sin(diptilt);  // sin(dipole tilt)
    tPS = sPS/cPS;  // tan(dipole tilt)
  
    posSM.x = posGSM->x*cPS - posGSM->z*sPS;  //  RBF EXPANSIONS ARE IN SM COORDINATES
    posSM.y = posGSM->y;                      //  ->  CONVERT X,Y,Z FROM GSW TO SM 
    posSM.z = posGSM->z*cPS + posGSM->x*sPS;
    
    Pdyn  = PARMOD[1];
    SymV  = PARMOD[2]/50.0;  // Sym-H divided by 50
    Xind  = PARMOD[3];
    ByIMF = PARMOD[4];
    FPD = sqrt(Pdyn/2.0)-1.0;
  
    BvecSM.x = 0.0; BvecSM.y = 0.0; BvecSM.z = 0.0;
    for (I=1; I<=1296; I++) {
      Parr.x = Info->XX[I];
      Parr.y = Info->YY[I];
      Parr.z = Info->ZZ[I];
      Marr.x = Parr.x;
      Marr.y = Parr.y;
      Marr.z = -Parr.z;
    
      DELTA_ZR = Info->RHBR[I]*tPS;
      DTHETA = -asin(DELTA_ZR)*Info->ST[I];
      sDT = sin(DTHETA);
      cDTM1 = cos(DTHETA)-1.0;
      DP.x = Parr.x*cDTM1+sDT*Info->ZCP[I];
      DP.y = Parr.y*cDTM1+sDT*Info->ZSP[I];
      DP.z = Parr.z*cDTM1-Info->RHO[I]*sDT;
      DM.x = Marr.x*cDTM1-sDT*Info->ZCP[I];
      DM.y = Marr.y*cDTM1-sDT*Info->ZSP[I];
      DM.z = Marr.z*cDTM1-Info->RHO[I]*sDT;
    
      CP = sqrt(pow(posSM.x-Parr.x-DP.x, 2) +
                pow(posSM.y-Parr.y-DP.y, 2) +
                pow(posSM.z-Parr.z-DP.z, 2) + D*D);    // RBF Ch_i+
      CM = sqrt(pow(posSM.x-Marr.x-DM.x, 2) +
                pow(posSM.y-Marr.y-DM.y, 2) +
                pow(posSM.z-Marr.z-DM.z, 2) + D*D);    // RBF Ch_i-

      Lgm_VecSub(&DCP, &posSM, &Parr);
      Lgm_VecSub(&DCP, &DCP, &DP);
      Lgm_VecDivideScalar(&DCP, CP);
      Lgm_VecSub(&DCM, &posSM, &Marr);
      Lgm_VecSub(&DCM, &DCM, &DM);
      Lgm_VecDivideScalar(&DCM, CM);
      Lgm_VecMult(&DCPsq, &DCP, &DCP);
      Lgm_VecDivideScalar(&DCPsq, CP);
      DCPsq.x = 1.0/CP - DCPsq.x;
      DCPsq.y = 1.0/CP - DCPsq.y;
      DCPsq.z = 1.0/CP - DCPsq.z;
      Lgm_VecMult(&DCMsq, &DCM, &DCM);
      Lgm_VecDivideScalar(&DCMsq, CM);
      DCMsq.x = 1.0/CM - DCMsq.x;
      DCMsq.y = 1.0/CM - DCMsq.y;
      DCMsq.z = 1.0/CM - DCMsq.z;
      DCPXY = -DCP.x*DCP.y/CP;
      DCPXZ = -DCP.x*DCP.z/CP;
      DCPYZ = -DCP.y*DCP.z/CP;
      DCMXY = -DCM.x*DCM.y/CM;
      DCMXZ = -DCM.x*DCM.z/CM;
      DCMYZ = -DCM.y*DCM.z/CM;
   
      TCP.x = posSM.z*DCP.y - posSM.y*DCP.z;
      TCP.y = posSM.x*DCP.z - posSM.z*DCP.x;
      TCP.z = posSM.y*DCP.x - posSM.x*DCP.y;
      TCM.x = posSM.z*DCM.y - posSM.y*DCM.z;
      TCM.y = posSM.x*DCM.z - posSM.z*DCM.x;
      TCM.z = posSM.y*DCM.x - posSM.x*DCM.y;
    
      PCP.x = 2.0*DCP.x-posSM.x*(DCPsq.y+DCPsq.z)+posSM.y*DCPXY+posSM.z*DCPXZ;
      PCP.y = 2.0*DCP.y-posSM.y*(DCPsq.x+DCPsq.z)+posSM.z*DCPYZ+posSM.x*DCPXY;
      PCP.z = 2.0*DCP.z-posSM.z*(DCPsq.x+DCPsq.y)+posSM.x*DCPXZ+posSM.y*DCPYZ;
      PCM.x = 2.0*DCM.x-posSM.x*(DCMsq.y+DCMsq.z)+posSM.y*DCMXY+posSM.z*DCMXZ;
      PCM.y = 2.0*DCM.y-posSM.y*(DCMsq.x+DCMsq.z)+posSM.z*DCMYZ+posSM.x*DCMXY;
      PCM.z = 2.0*DCM.z-posSM.z*(DCMsq.x+DCMsq.y)+posSM.x*DCMXZ+posSM.y*DCMYZ;

      Lgm_VecAdd(&CTarr, &TCP, &TCM); 
      Lgm_VecMultiplyScalar(&CTarr, cPS);
      Lgm_VecSub(&STarr, &TCP, &TCM); 
      Lgm_VecMultiplyScalar(&STarr, sPS);
      Lgm_VecSub(&CParr, &PCP, &PCM); 
      Lgm_VecMultiplyScalar(&CParr, cPS);
      Lgm_VecAdd(&SParr, &PCP, &PCM); 
      Lgm_VecMultiplyScalar(&SParr, sPS);
    
    // -----------------   TOTAL FIELD:    -----------------------------------
      ACT = Info->A[I]+Info->A[I+5184]*FPD+Info->A[I+10368]*SymV+Info->A[I+15552]*Xind;
      AST = Info->A[I+1296]+Info->A[I+6480]*FPD+Info->A[I+11664]*SymV+Info->A[I+16848]*Xind;
      AT  = Info->A[I+20736]*ByIMF;
      ACP = Info->A[I+2592]+Info->A[I+7776]*FPD+Info->A[I+12960]*SymV+Info->A[I+18144]*Xind;
      ASP = Info->A[I+3888]+Info->A[I+9072]*FPD+Info->A[I+14256]*SymV+Info->A[I+19440]*Xind;
      AP  = Info->A[I+22032]*ByIMF;
      BvecSM.x = BvecSM.x + CTarr.x*ACT + STarr.x*AST + (TCP.x-TCM.x)*AT + CParr.x*ACP + SParr.x*ASP + (PCP.x+PCM.x)*AP;
      BvecSM.y = BvecSM.y + CTarr.y*ACT + STarr.y*AST + (TCP.y-TCM.y)*AT + CParr.y*ACP + SParr.y*ASP + (PCP.y+PCM.y)*AP;
      BvecSM.z = BvecSM.z + CTarr.z*ACT + STarr.z*AST + (TCP.z-TCM.z)*AT + CParr.z*ACP + SParr.z*ASP + (PCP.z+PCM.z)*AP;
    }
    Lgm_Convert_Coords(&BvecSM, BvecGSM, SM_TO_GSM, ctrans);
    return(1);
}

int Lgm_B_TA16( Lgm_Vector *rGSM, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {
    // Evaluate Tsyganenko & Andreeva 2016 model
    Lgm_Vector      Bext, Bint;
    double          PARMOD[11];
    int             retval;

    PARMOD[1] = Info->TA16_Info.Pdyn;
    PARMOD[2] = Info->TA16_Info.SymHc_avg;
    PARMOD[3] = Info->TA16_Info.Xind_avg;
    PARMOD[4] = Info->TA16_Info.By_avg;
    retval = TA2016(rGSM, PARMOD, Info->c, &Bext, &Info->TA16_Info);
    switch ( Info->InternalModel ){
        case LGM_CDIP:
                        Lgm_B_cdip(rGSM, &Bint, Info);
                        break;
        case LGM_EDIP:
                        Lgm_B_edip(rGSM, &Bint, Info);
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf(rGSM, &Bint, Info);
                        break;
        default:
                        fprintf(stderr, "Lgm_B_TA16: Unknown internal model (%d)\n", Info->InternalModel );
                        break;
    }

    // Add requested internal field model to obtain total field
    Lgm_VecAdd(B, &Bint, &Bext);

    ++Info->nFunc;  // increment counter of function calls

    return(1);

}
