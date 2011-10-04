#include "Lgm/Lgm_CTrans.h"

/* This is where I store the date */
Lgm_DateTime AEdatetime[4096] ;
int AE_index[4096], AL_index[4096], AU_index[4096], AO_index[4096], index_length ;
char AEfilename[1024] ;

#ifdef __cplusplus
extern "C" {
#endif

/* Function definitions */
double Lgm_get_AE(int,int,int,int,int,int) ;
double Lgm_get_AL(int,int,int,int,int,int) ;
double Lgm_get_AU(int,int,int,int,int,int) ;
double Lgm_get_AO(int,int,int,int,int,int) ;
double Lgm_get_index(int,int,int,int,int,int,int*) ;
void Lgm_Read_AEfile(int,int,int) ;

#ifdef __cplusplus
}
#endif
