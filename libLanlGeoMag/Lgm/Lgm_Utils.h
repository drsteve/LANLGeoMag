#ifndef LGM_UTILS_H
#define LGM_UTILS_H

void    Lgm_LogSpace( double start, double stop, int num, double *array );
void    Lgm_LinSpace(double start, double stop, int num, double *array);
int     Lgm_Bisect(double *data, double value, int len);




#endif
