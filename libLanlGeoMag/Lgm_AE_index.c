#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_statistics_double.h>
#include "Lgm/Lgm_AE_index.h"

/**
 * \brief
 * Return the AE index at a given year, month, day, hour, minute and sec.
 *
 * Data is stored in /n/projects/lanl/geo/Data/AE from 1975 to 2011
 *
 * \details
 * 
 * \param[in] year (1975-2011)
 * \param[in] month (1-12)
 * \param[in] day (1-31)
 * \param[in] hour (1-24)
 * \param[in] minute (1-60)
 * \param[in] second (1-60)
 *
 * \return AE index; other functions return AL, AU, AO
 *
 * \author  C. Jeffery
 * \date 2011
 *
 *
 */

double Lgm_get_AE(int year, int month, int day, int hour, int min, int sec) {

  int* index_ptr ;
  double AE ;

  index_ptr = AE_index ;

  AE = Lgm_get_index(year,month,day,hour,min,sec,index_ptr) ;

  return(AE) ;
}

/*
 * \brief
 * Similar to Lgm_get_AE (above) but return AL index
 */
double Lgm_get_AL(int year, int month, int day, int hour, int min, int sec) {

  int* index_ptr ;
  double AL ;

  index_ptr = AL_index ;

  AL = Lgm_get_index(year,month,day,hour,min,sec,index_ptr) ;

  return(AL) ;
}

/*
 * \brief
 * Similar to Lgm_get_AE (above) but return AU index
 */
double Lgm_get_AU(int year, int month, int day, int hour, int min, int sec) {

  int* index_ptr ;
  double AUo ;

  index_ptr = AU_index ;

  AUo = Lgm_get_index(year,month,day,hour,min,sec,index_ptr) ;

  return(AUo) ;
}

/*
 * \brief
 * Similar to Lgm_get_AE (above) but return AO index
 */
double Lgm_get_AO(int year, int month, int day, int hour, int min, int sec) {

  int* index_ptr ;
  double AO ;

  index_ptr = AO_index ;

  AO = Lgm_get_index(year,month,day,hour,min,sec,index_ptr) ;

  return(AO) ;
}

/*
 * \brief
 * Generic index retrieval
 */
double Lgm_get_index(int year, int month, int day, int hour, int min, int sec, int* index_ptr) {

  int cntr,min_index,next_index ;
  double fYear,t1,t2,indexfit ;
  double* fYear_array, *Data_array ;
  double      MJD, Prev_MJD, Next_MJD, Prev_UT, Next_UT;
  long int    total_length=0, Prev_Date, Next_Date;
  int         Prev_Year, Prev_Month, Prev_Day, Next_Year, Next_Month, Next_Day;  
  Lgm_CTrans  *c = Lgm_init_ctrans(0);
  

  fYear_array = malloc(4500*sizeof(double)) ;
  Data_array = malloc(4500*sizeof(double)) ;
  
  MJD = Lgm_MJD( year, month, day, 12.0, LGM_TIME_SYS_UTC, c );

  Prev_MJD = MJD-1.0;
  Lgm_mjd_to_ymdh( Prev_MJD, &Prev_Date, &Prev_Year, &Prev_Month, &Prev_Day, &Prev_UT );
  //printf("Previous day with year=%d, month=%d, day=%d\n",Prev_Year, Prev_Month, Prev_Day);

  Next_MJD = MJD+1.0;
  Lgm_mjd_to_ymdh( Next_MJD, &Next_Date, &Next_Year, &Next_Month, &Next_Day, &Next_UT );
  //printf("Next day with year=%d, month=%d, day=%d\n",Next_Year, Next_Month, Next_Day);
  
  //Read in previous day
  Lgm_Read_AEfile(Prev_Year, Prev_Month, Prev_Day) ;
  for(cntr=0; cntr<index_length; cntr++)
    {
      fYear_array[cntr] = AEdatetime[cntr].fYear ;
      Data_array[cntr] = index_ptr[cntr];
    }
  total_length+=index_length;
  
  //Read in current day
  Lgm_Read_AEfile(year, month, day) ;
  for(cntr=0; cntr<index_length; cntr++)
    {
      fYear_array[cntr+total_length] = AEdatetime[cntr].fYear ;
      Data_array[cntr+total_length] = index_ptr[cntr];
    }
  total_length+=index_length;

  //Read in next day
  Lgm_Read_AEfile(Next_Year, Next_Month, Next_Day) ;
  for(cntr=0; cntr<index_length; cntr++)
    {
      fYear_array[cntr+total_length] = AEdatetime[cntr].fYear ;
      Data_array[cntr+total_length] = index_ptr[cntr];
    }
  total_length+=index_length  ;  

  fYear = year + month/12. + day/(12.*31.) + hour/(12.*31.*24.) + min/(12.*31.*24.*60.) + 
    sec/(12.*31.*24.*60.*60.) ;

  for(cntr=0; cntr<total_length; cntr++)
    {
      fYear_array[cntr] = fabs(fYear_array[cntr] - fYear) ;
    }

  /* figure out my best match */

  min_index = gsl_stats_min_index(fYear_array,1,total_length) ;
  
  /* do a simple linear interpolation */

  t1 = fYear_array[min_index] ;

  if(fYear_array[min_index-1]>fYear_array[min_index+1])
    next_index = min_index+1 ;
  else
    next_index = min_index-1 ;

  t2 = fYear_array[next_index] ;

  /* simple weighting function */

  indexfit = (t2*Data_array[min_index] + t1*Data_array[next_index])/(t2+t1) ;

  free(fYear_array) ;

  return(indexfit) ;
}

/*
 * \brief
 * Create file name
 */
void Lgm_get_filename(char* year, char* month, char* day) {

  int month_, day_ ;

  month_ = strtol(month,&month,0) ;
  day_ = strtol(day,&day,0) ;
  
  sprintf(AEfilename,"/n/projects/lanl/geo/Data/AE/%4.4s/AE_%4.4s%2.2d%2.2d.dat", year, year, month_, day_) ;

}
  
/*
 * \brief
 * Read in AE index
 */
void Lgm_Read_AEfile(int year, int month, int day) {

  FILE *fp ;
  char buffer[1024], UTC[1024] ;
  int cntr,hold_second ;

  sprintf(AEfilename,"/n/projects/lanl/geo/Data/AE/%4.4d/AE_%4.4d%2.2d%2.2d.dat", year, year, month, day) ;

  /* Open file */

  fp = fopen(AEfilename, "r") ;
  
  if(fp == NULL) {

    printf("Can't open file=%s\n",AEfilename) ;

    if(year<1975)
      printf("Year must be between 1975 and 2011\n") ;
    else if(year>2011)
      printf("Need to update /n/projects/lanl/geo/Data/AE with new data\n") ;
    else if(month<1)
      printf("Month must be between 1 and 12\n") ;
    else if(month>12)
      printf("Month must be between 1 and 12\n") ;
    else if(day<1)
      printf("Day must be between 1 and 31\n") ;
    else if(day>31)
      printf("Day must be between 1 and 31\n") ;
    else
      printf("Data file appears to be missing\n") ;

    exit(1) ;
  }

  /* advance to record data */

  do {
    fgets(buffer,1024,fp) ;
  }
  while(buffer[0]=='#') ;

  cntr = 0 ;

  do {
    sscanf(buffer,"%s %d %d %d %d\n", &UTC, &AE_index[cntr], &AL_index[cntr], &AU_index[cntr], 
	   &AO_index[cntr]) ;

    sscanf(UTC,"%d-%d-%dT%d:%d:%d", &AEdatetime[cntr].Year, &AEdatetime[cntr].Month, 
	   &AEdatetime[cntr].Day, &AEdatetime[cntr].Hour, &AEdatetime[cntr].Minute, 
	   &hold_second) ;

    /* convert second into double */

    AEdatetime[cntr].Second = hold_second ;

    AEdatetime[cntr].fYear = AEdatetime[cntr].Year + AEdatetime[cntr].Month/12. + 
      AEdatetime[cntr].Day/(12.*31.) + AEdatetime[cntr].Hour/(12.*31.*24.) + 
      AEdatetime[cntr].Minute/(12.*31.*24.*60.) + AEdatetime[cntr].Second/(12.*31.*24.*60.*60.) ;

    cntr += 1 ;

  }
  while (fgets(buffer,1024,fp) != NULL) ;

  index_length = cntr ;

  fclose(fp) ;
}
