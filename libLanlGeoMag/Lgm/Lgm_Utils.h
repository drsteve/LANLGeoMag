#ifndef LGM_UTILS_H
#define LGM_UTILS_H

void    Lgm_LogSpace( double start, double stop, long num, double *array );
void    Lgm_LinSpace(double start, double stop, long num, double *array);
long     Lgm_Bisect(double *data, double value, long len);




#endif
