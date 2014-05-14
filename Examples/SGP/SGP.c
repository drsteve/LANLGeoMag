



#include <stdio.h>
#include <stdlib.h>

#include <Lgm_Sgp.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagModelInfo.h>


int main(void) {
  char            *filename ="iss.tle";
  int             nTle=0; // this is not the number of TLE but the index in the tle array
  _SgpTLE         *tle;  // pointer to a struct
  Lgm_CTrans      *c;
  // create a place to hold an array of SGP trackers (1 element)
  _SgpInfo        *s;
  double           Lat, Lon, r, minutes, tsince=0, JD, tUT;
  int              tYear, tMonth, tDay;
  long int         tDate, i;
  Lgm_Vector       Ugei;

  // create a place to hold an array of TLEs (1 element)
  tle = (_SgpTLE *)calloc( 1, sizeof(_SgpTLE) );

  if (!LgmSgp_ReadTlesFromFile( filename, &nTle, tle, 4)){
    printf("TLE not parsed!\n");
    free(tle);
    return(-1);
  } else {
    printf("TLE read and parsed.\n");
  }
  
  s = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
  c = Lgm_init_ctrans( 0 );

  printf("Line0: %s\n", tle[0].Line0);
  printf("Line1: %s\n", tle[0].Line1);
  printf("Line2: %s\n", tle[0].Line2);
  printf("\n");

    /*
     * All the TLEs have their own epoch times0 in them. And the propagator (sgp4)
     * uses the "time since (in minutes)". So for a given time of interest, we need to
     * compute the tsince needed.
     */
  // JD=2454771.014400;
  JD = Lgm_JD( 2008, 10, 31, 12.3456, LGM_TIME_SYS_UTC, c );
  minutes =  12.3456;
  printf("JD( 2008, 10, 31, 12.3456, c  ) = %lf\n",  JD);

  Lgm_jd_to_ymdh( JD, &tDate, &tYear, &tMonth, &tDay, &tUT );
  Lgm_Set_Coord_Transforms( tDate, tUT, c );

  tsince = (JD - tle[0].JD)*1440.0;
  printf("tsince=%lf (minutes from TLE epoch)\n", tsince);
  printf("\n");
  printf("Compute the ground track for 2 hours:\n");
  // that is 120 minutes / 5 minutes per point
  
  // initialize the propagator
  LgmSgp_SGP4_Init( s, tle );
  for (i=0; i<120/3; i++) {
    JD = Lgm_JD( 2008, 10, 31, minutes, LGM_TIME_SYS_UTC, c );
    Lgm_jd_to_ymdh( JD, &tDate, &tYear, &tMonth, &tDay, &tUT );
    tsince = (JD - tle[0].JD)*1440.0;
    // coord transform as time depend to much be resetup
    Lgm_Set_Coord_Transforms( tDate, tUT, c );
    // do the propagation
    LgmSgp_SGP4( tsince, s );
    // make into Re numbers
    Ugei.x = s->X/Re; 
    Ugei.y = s->Y/Re; 
    Ugei.z = s->Z/Re;
    // get lat and lon
    Lgm_CartToSphCoords(&Ugei, &Lat, &Lon, &r);
    printf("\t%ld %lf\tLat:%lf\tLon:%lf\n", tDate, tUT, Lat, Lon);
    
    minutes += 3./60.; // add three minutes
  }

  Lgm_free_ctrans( c ); // free the structure
  free(s);
  free(tle);
  return(0);
}


