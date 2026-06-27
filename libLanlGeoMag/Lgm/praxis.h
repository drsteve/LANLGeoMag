#ifndef LGM_PRAXIS_H
#define LGM_PRAXIS_H

double *allocate_real_vector(int l, int u);
double **allocate_real_matrix(int lr, int ur, int lc, int uc);
void free_real_vector(double* v, int l);
void free_real_matrix(double** m, int lr, int ur, int lc);
void inivec(int l, int u, double* a, double x);
void inimat(int lr, int ur, int lc, int uc, double** a, double x);
void dupvec(int l, int u, int shift, double* a,  double* b);
void dupcolvec(int l, int u, int j, double** a, double* b);
void dupmat(int l, int u, int i, int j, double** a, double** b);
double mattam(int l, int u, int i, int j, double** a, double** b);
double matmat(int l, int u, int i, int j, double** a, double** b);
void hshreabid(double** a, int m, int n, double* d, double* b, double* em);
void mulcol(int l, int u, int i, int j, double** a, double** b, double x);
void mulrow(int l, int u, int i, int j, double** a, double** b, double x);
double tammat(int l, int u, int i, int j, double** a, double** b);
double vecvec(int l, int u, int shift, double* a, double* b);
void elmrow(int l, int u, int i, int j, double** a, double** b, double x);
void elmcol(int l, int u, int i, int j, double** a, double** b, double x);
void ichrowcol(int l, int u, int i, int j, double** a);
void elmveccol(int l, int u, int i, double* a, double** b, double x);
int qrisngvaldec(double** a, int m, int n, double* val, double** v, double* em);
int qrisngvaldecbid(double* d, double* b, int m, int n, double** u, double** v, double* em);
void psttfmmat(double** a, int n, double** v, double* b);
void pretfmmat(double **a, int m, int n, double* d);
void praxismin(int j, int nits, double* d2, double* x1, double* f1, int fk,
  int n, double* x, double** v, double* qa, double* qb, double* qc,
  double qd0, double qd1, double* q0, double* q1, int* nf, int* nl,
  double* fx, double m2, double m4, double dmin, double ldt, double reltol,
  double abstol, double small, double h, double (*funct)(double*, void *data),
  int* data);
void rotcol(int l, int u, int i, int j, double** a, double c, double s);
double praxisflin(double l, int j, int n, double* x, double** v,
    double *qa, double* qb, double* qc, double qd0,double qd1,
    double* q0, double* q1, int* nf, double(*f)(double*, void*), int* data);
void praxis( int n, double *x, int *data, double (*funct)(double *, void *data),
    double *in, double *out);

#endif
