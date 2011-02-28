/*
 *
 *    praxis.c -- not sure where the original source of this code is.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <values.h>  // removed by BAL 28Feb2011, my mac doesn't have this header is it really needed?
#include <float.h>

#ifdef LINUX
#include <values.h>
#endif

#ifdef CYGWIN
#define DBL_MIN (1e-37)
#endif

#ifdef MSWIN
#include <limits.h>
#include <float.h>
#endif




/* 
 *  Some defines that may not be known
 *  by all gcc compilers...
 */
#ifndef RAND_MAX
  /* 2^15 - 1 */
#define RAND_MAX (32767.0)
#endif




	double *allocate_real_vector();
	double **allocate_real_matrix();
	void free_real_vector();
	void free_real_matrix();
	void inivec();
	void inimat();
	void dupvec();
	void dupmat();
	void dupcolvec();
	void mulrow();
	void mulcol();
	double vecvec();
	double tammat();
	double mattam();
	void ichrowcol();
	void elmveccol();
	int qrisngvaldec();
	void praxismin();



void praxis( int n, double *x, int *data, double (*funct)(double *, int *data), double *in, double *out) {

	int illc,i,j,k,k2,nl,maxf,nf,kl,kt,ktm,emergency;
	double s,sl,dn,dmin,fx,f1,lds,ldt,sf,df,qf1,qd0,qd1,qa,qb,qc,m2,m4,
			small,vsmall,large,vlarge,scbd,ldfac,t2,macheps,reltol,
			abstol,h,**v,*d,*y,*z,*q0,*q1,**a,em[8],l;

	/*
	 *  Seed random number generator
	 */
#ifdef MSWIN
	srand(34084320);
#else
	srand48(34084320);
#endif
	
//	for (i=0; i<8; ++i) x[i+1] = (double)data->x[i];
	d=allocate_real_vector(1,n);
	y=allocate_real_vector(1,n);
	z=allocate_real_vector(1,n);
	q0=allocate_real_vector(1,n);
	q1=allocate_real_vector(1,n);
	v=allocate_real_matrix(1,n,1,n);
	a=allocate_real_matrix(1,n,1,n);

    //  heuristic numbers:
    //
    //  If the axes may be badly scaled (which is to be avoided if
    //  possible), then set scbd = 10.  otherwise set scbd=1.
    //
    //  If the problem is known to be ill-conditioned, set ILLC = true.
    //
    //  KTM is the number of iterations without improvement before the
    //  algorithm terminates.  KTM = 4 is very cautious; usually KTM = 1
    //  is satisfactory.
    //

	macheps=in[0];
	reltol=in[1];
	abstol=in[2];
	maxf=in[5];
	h=in[6];
	scbd=in[7];
	ktm=in[8];
	illc = in[9] < 0.0;
	small=macheps*macheps;
	vsmall=small*small;
	large=1.0/small;
	vlarge=1.0/vsmall;
	m2=reltol;
	m4=sqrt(m2);
	srand(1);
	ldfac = (illc ? 0.1 : 0.01);
	kt=nl=0;
	nf=1;
	out[3]=qf1=fx=(*funct)(x, data);
	abstol=t2=small+fabs(abstol);
	dmin=small;
	if (h < abstol*100.0) h=abstol*100;
	ldt=h;
	inimat(1,n,1,n,v,0.0);
	for (i=1; i<=n; i++) v[i][i]=1.0;
	d[1]=qd0=qd1=0.0;
	dupvec(1,n,0,q1,x);
	inivec(1,n,q0,0.0);
	emergency=0;

	while (1) {
		sf=d[1];
		d[1]=s=0.0;
		praxismin(1,2,&(d[1]),&s,&fx,0,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct, data);
		if (s <= 0.0) mulcol(1,n,1,1,v,v,-1.0);
		if (sf <= 0.9*d[1] || 0.9*sf >= d[1]) inivec(2,n,d,0.0);
		for (k=2; k<=n; k++) {
			dupvec(1,n,0,y,x);
			sf=fx;
			illc = (illc || kt > 0);
			while (1) {
				kl=k;
				df=0.0;
				if (illc) {
					/* random stop to get off resulting valley */
					for (i=1; i<=n; i++) {
						s=z[i]=(0.1*ldt+t2*pow(10.0,kt))*
#ifdef MSWIN
									((double)(rand())/RAND_MAX-0.5);
#else
									(drand48()-0.5);
#endif
						elmveccol(1,n,i,x,v,s);
					}
					fx=(*funct)(x, data);
					nf++;
				}
				for (k2=k; k2<=n; k2++) {
					sl=fx;
					s=0.0;
					praxismin(k2,2,&(d[k2]),&s,&fx,0,
						n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
						&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct, data);
					s = illc ? d[k2]*(s+z[k2])*(s+z[k2]) : sl-fx;
					if (df < s) {
						df=s;
						kl=k2;
					}
				}
				if (!illc && df < fabs(100.0*macheps*fx))
					illc=1;
				else
					break;
			}
			for (k2=1; k2<=k-1; k2++) {
				s=0.0;
				praxismin(k2,2,&(d[k2]),&s,&fx,0,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct, data);
			}
			f1=fx;
			fx=sf;
			lds=0.0;
			for (i=1; i<=n; i++) {
				sl=x[i];
				x[i]=y[i];
				y[i] = sl -= y[i];
				lds += sl*sl;
			}
			lds=sqrt(lds);
			if (lds > small) {
				for (i=kl-1; i>=k; i--) {
					for (j=1; j<=n; j++) v[j][i+1]=v[j][i];
					d[i+1]=d[i];
				}
				d[k]=0.0;
				dupcolvec(1,n,k,v,y);
				mulcol(1,n,k,k,v,v,1.0/lds);
				praxismin(k,4,&(d[k]),&lds,&f1,1,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct, data);
				if (lds <= 0.0) {
					lds = -lds;
					mulcol(1,n,k,k,v,v,-1.0);
				}
			}
			ldt *= ldfac;
			if (ldt < lds) ldt=lds;
			t2=m2*sqrt(vecvec(1,n,0,x,x))+abstol;
			kt = (ldt > 0.5*t2) ? 0 : kt+1;
			if (kt > ktm) {
				out[1]=0.0;
				emergency=1;
			}
		}
		if (emergency) break;
		/* quad */
		s=fx;
		fx=qf1;
		qf1=s;
		qd1=0.0;
		for (i=1; i<=n; i++) {
			s=x[i];
			x[i]=l=q1[i];
			q1[i]=s;
			qd1 += (s-l)*(s-l);
		}
		l=qd1=sqrt(qd1);
		s=0.0;
		if ((qd0*qd1 > DBL_MIN) && (nl >=3*n*n)) {
			praxismin(0,2,&s,&l,&qf1,1,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct, data);
			qa=l*(l-qd1)/(qd0*(qd0+qd1));
			qb=(l+qd0)*(qd1-l)/(qd0*qd1);
			qc=l*(l+qd0)/(qd1*(qd0+qd1));
		} else {
			fx=qf1;
			qa=qb=0.0;
			qc=1.0;
		}
		qd0=qd1;
		for (i=1; i<=n; i++) {
			s=q0[i];
			q0[i]=x[i];
			x[i]=qa*s+qb*x[i]+qc*q1[i];
		}
		/* end of quad */
		dn=0.0;
		for (i=1; i<=n; i++) {
			d[i]=1.0/sqrt(d[i]);
			if (dn < d[i]) dn=d[i];
		}
		for (j=1; j<=n; j++) {
			s=d[j]/dn;
			mulcol(1,n,j,j,v,v,s);
		}
		if (scbd > 1.0) {
			s=vlarge;
			for (i=1; i<=n; i++) {
				sl=z[i]=sqrt(mattam(1,n,i,i,v,v));
				if (sl < m4) z[i]=m4;
				if (s > sl) s=sl;
			}
			for (i=1; i<=n; i++) {
				sl=s/z[i];
				z[i]=1.0/sl;
				if (z[i] > scbd) {
					sl=1.0/scbd;
					z[i]=scbd;
				}
				mulrow(1,n,i,i,v,v,sl);
			}
		}
		for (i=1; i<=n; i++) ichrowcol(i+1,n,i,i,v);
		em[0]=em[2]=macheps;
		em[4]=10*n;
		em[6]=vsmall;
		dupmat(1,n,1,n,a,v);
		if (qrisngvaldec(a,n,n,d,v,em) != 0) {
			out[1]=2.0;
			emergency=1;
		}
		if (emergency) break;
		if (scbd > 1.0) {
			for (i=1; i<=n; i++) mulrow(1,n,i,i,v,v,z[i]);
			for (i=1; i<=n; i++) {
				s=sqrt(tammat(1,n,i,i,v,v));
				d[i] *= s;
				s=1.0/s;
				mulcol(1,n,i,i,v,v,s);
			}
		}
		for (i=1; i<=n; i++) {
			s=dn*d[i];
			d[i] = (s > large) ? vsmall :
						((s < small) ? vlarge : 1.0/(s*s));
		}
		/* sort */
		for (i=1; i<=n-1; i++) {
			k=i;
			s=d[i];
			for (j=i+1; j<=n; j++)
				if (d[j] > s) {
					k=j;
					s=d[j];
				}
			if (k > i) {
				d[k]=d[i];
				d[i]=s;
				for (j=1; j<=n; j++) {
					s=v[j][i];
					v[j][i]=v[j][k];
					v[j][k]=s;
				}
			}
		}
		/* end of sort */
		dmin=d[n];
		if (dmin < small) dmin=small;
		illc = (m2*d[1]) > dmin;
		if (nf >= maxf) {
			out[1]=1.0;
			break;
		}
	}
	out[2]=fx;
	out[4]=nf;
	out[5]=nl;
	out[6]=ldt;
	free_real_vector(d,1);
	free_real_vector(y,1);
	free_real_vector(z,1);
	free_real_vector(q0,1);
	free_real_vector(q1,1);
	free_real_matrix(v,1,n,1);
	free_real_matrix(a,1,n,1);

//	for (i=0; i<40; ++i) data->x[i] = (double)x[i+1];

}

void praxismin(j, nits, d2, x1, f1, fk, n, x, v, qa, qb, qc, qd0, qd1, q0, q1, nf, nl, 
			fx, m2, m4, dmin, ldt, reltol, abstol, small, h, funct, data)
int 	j;
int 	nits;
double 	*d2;
double 	*x1;
double 	*f1;
int 	fk;
int 	n;
double 	x[];
double 	**v;
double 	*qa;
double 	*qb;
double 	*qc;
double 	qd0;
double 	qd1;
double 	q0[];
double 	q1[];
int 	*nf;
int 	*nl;
double 	*fx;
double 	m2;
double	m4;
double 	dmin;
double 	ldt;
double 	reltol;
double 	abstol;
double 	small;
double 	h;
double   (*funct)(double *, int *data);
int *data;



{
	/* this function is internally used by PRAXIS */

	double praxisflin();
	int k,dz,loop;
	double x2,xm,f0,f2,fm,d1,t2,s,sf1,sx1;

	sf1 = *f1;
	sx1 = *x1;
	k=0;
	xm=0.0;
	f0 = fm = *fx;
	dz = *d2 < reltol;
	s=sqrt(vecvec(1,n,0,x,x));
	t2=m4*sqrt(fabs(*fx)/(dz ? dmin : *d2)+s*ldt)+m2*ldt;
	s=s*m4+abstol;
	if (dz && (t2 > s)) t2=s;
	if (t2 < small) t2=small;
	if (t2 > 0.01*h) t2=0.01*h;
	if (fk && (*f1 <= fm)) {
		xm = *x1;
		fm = *f1;
	}
	if (!fk || (fabs(*x1) < t2)) {
		*x1 = (*x1 > 0.0) ? t2 : -t2;
		*f1=praxisflin(*x1,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct, data);
	}
	if (*f1 <= fm) {
		xm = *x1;
		fm = *f1;
	}
	loop=1;
	while (loop) {
		if (dz) {
			/* evaluate praxisflin at another point and
				estimate the second derivative */
			x2 = (f0 < *f1) ? -(*x1) : (*x1)*2.0;
			f2=praxisflin(x2,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct, data);
			if (f2 <= fm) {
				xm=x2;
				fm=f2;
			}
			*d2=(x2*((*f1)-f0)-(*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
		}
		/* estimate first derivative at 0 */
		d1=((*f1)-f0)/(*x1)-(*x1)*(*d2);
		dz=1;
		x2 = (*d2 <= small) ? ((d1 < 0.0) ? h : -h) : -0.5*d1/(*d2);
		if (fabs(x2) > h)	x2 = (x2 > 0.0) ? h : -h;
		while (1) {
			f2=praxisflin(x2,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct, data);
			if (k < nits && f2 > f0) {
				k++;
				if (f0 < *f1 && (*x1)*x2 > 0.0) break;
				x2=0.5*x2;
			} else {
				loop=0;
				break;
			}
		}
	}
	(*nl)++;
	if (f2 > fm)
		x2=xm;
	else
		fm=f2;
	*d2 = (fabs(x2*(x2-(*x1))) > small) ?
				((x2*((*f1)-f0)-(*x1)*(fm-f0))/((*x1)*x2*((*x1)-x2))) :
				((k > 0) ? 0.0 : *d2);
	if (*d2 <= small) *d2=small;
	*x1=x2;
	*fx=fm;
	if (sf1 < *fx) {
		*fx=sf1;
		*x1=sx1;
	}
	if (j > 0) elmveccol(1,n,j,x,v,*x1);
}

double praxisflin(l, j, n, x, v, qa, qb, qc, qd0, qd1, q0, q1, nf, funct, data)

double 	l; 
int 	j; 
int 	n; 
double 	x[]; 
double 	**v;
double 	*qa; 
double 	*qb; 
double 	*qc; 
double 	qd0; 
double 	qd1; 
double 	q0[]; 
double 	q1[]; 
int 	*nf; 
double  (*funct)(double *, int *);
int     *data;

{
	/* this function is internally used by PRAXISMIN */

	int i;
	double *t,result;

	t=allocate_real_vector(1,n);
	if (j > 0)
		for (i=1; i<=n; i++) t[i]=x[i]+l*v[i][j];
	else {
		/* search along parabolic space curve */
		*qa=l*(l-qd1)/(qd0*(qd0+qd1));
		*qb=(l+qd0)*(qd1-l)/(qd0*qd1);
		*qc=l*(l+qd0)/(qd1*(qd0+qd1));
		for (i=1; i<=n; i++) t[i]=(*qa)*q0[i]+(*qb)*x[i]+(*qc)*q1[i];
	}
	(*nf)++;
	result=(*funct)(t, data);
	free_real_vector(t,1);
	return result;
}

void dupcolvec(l, u, j, a, b)
int 	l;
int 	u;
int 	j;
double 	**a;
double 	b[];
{
	for (; l<=u; l++) a[l][j]=b[l];
}

void dupmat(l, u, i, j, a, b)
int 	l;
int 	u;
int 	i;
int 	j;
double 	**a;
double 	**b;
{
	int k;

	for (; l<=u; l++)
		for (k=i; k<=j; k++) a[l][k]=b[l][k];
}



void dupvec(l, u, shift, a, b)
int l; 
int u; 
int shift; 
double a[]; 
double b[];

{
	for (; l<=u; l++) a[l]=b[l+shift];
}



void elmveccol(l, u, i, a, b, x)
int l; 
int u; 
int i; 
double a[]; 
double **b; 
double x;
{
	for (; l<=u; l++) a[l] += b[l][i]*x;
}



void ichrowcol(l, u, i, j, a)
int l; 
int u; 
int i; 
int j; 
double **a;
{
	double r;

	for (; l<=u; l++) {
		r=a[i][l];
		a[i][l]=a[l][j];
		a[l][j]=r;
	}
}



void inimat(lr, ur, lc, uc, a, x)
int lr;
int ur;
int lc;
int uc;
double **a;
double x;
{
	int j;

	for (; lr<=ur; lr++)
		for (j=lc; j<=uc; j++) a[lr][j]=x;
}


void inivec(l, u, a, x)
int l;
int u;
double a[];
double x;
{
	for (; l<=u; l++) a[l]=x;
}




double mattam(l, u, i, j, a, b)
int l;
int u;
int i;
int j;
double **a;
double **b;
{
	int k;
	double s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[i][k]*b[j][k];
	return (s);
}



void mulcol(l, u, i, j, a, b, x)
int l;
int u;
int i;
int j;
double **a;
double **b;
double x;
{
	for (; l<=u; l++) a[l][i]=b[l][j]*x;
}


void mulrow(l, u, i, j, a, b, x)
int l;
int u;
int i;
int j;
double **a;
double **b;
double x;
{
	for (; l<=u; l++) a[i][l]=b[j][l]*x;
}




double tammat(l, u, i, j, a, b)
int l;
int u;
int i;
int j;
double **a;
double **b;
{
	int k;
	double s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[k][i]*b[k][j];
	return (s);
}



double vecvec(l, u, shift, a, b)
int l;
int u;
int shift;
double a[];
double b[];
{
	int k;
	double s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[k]*b[k+shift];
	return (s);
}



int qrisngvaldec(a, m, n, val, v, em)
double **a;
int m;
int n;
double val[];
double **v;
double em[];
{
	double *allocate_real_vector();
	void free_real_vector();
	void hshreabid();
	void psttfmmat();
	void pretfmmat();
	int qrisngvaldecbid();
	int i;
	double *b;

	b=allocate_real_vector(1,n);
	hshreabid(a,m,n,val,b,em);
	psttfmmat(a,n,v,b);
	pretfmmat(a,m,n,val);
	i=qrisngvaldecbid(val,b,m,n,a,v,em);
	free_real_vector(b,1);
	return i;
}


int qrisngvaldecbid(d, b, m, n, u, v, em)
double d[];
double b[];
int m;
int n;
double **u;
double **v;
double em[];
{
	void rotcol();
	int n0,n1,k,k1,i,i1,count,max,rnk;
	double tol,bmax,z,x,y,g,h,f,c,s,min;

	tol=em[2]*em[1];
	count=0;
	bmax=0.0;
	max=em[4];
	min=em[6];
	rnk=n0=n;
	do {
		k=n;
		n1=n-1;
		while (1) {
			k--;
			if (k <= 0) break;
			if (fabs(b[k]) >= tol) {
				if (fabs(d[k]) < tol) {
					c=0.0;
					s=1.0;
					for (i=k; i<=n1; i++) {
						f=s*b[i];
						b[i] *= c;
						i1=i+1;
						if (fabs(f) < tol) break;
						g=d[i1];
						d[i1]=h=sqrt(f*f+g*g);
						c=g/h;
						s = -f/h;
						rotcol(1,m,k,i1,u,c,s);
					}
					break;
				}
			} else {
				if (fabs(b[k]) > bmax) bmax=fabs(b[k]);
				break;
			}
		}
		if (k == n1) {
			if (d[n] < 0.0) {
				d[n] = -d[n];
				for (i=1; i<=n0; i++) v[i][n] = -v[i][n];
			}
			if (d[n] <= min) rnk--;
			n=n1;
		} else {
			count++;
			if (count > max) break;
			k1=k+1;
			z=d[n];
			x=d[k1];
			y=d[n1];
			g = (n1 == 1) ? 0.0 : b[n1-1];
			h=b[n1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=sqrt(f*f+1.0);
			f=((x-z)*(x+z)+h*(y/((f < 0.0) ? f-g : f+g)-h))/x;
			c=s=1.0;
			for (i=k1+1; i<=n; i++) {
				i1=i-1;
				g=b[i1];
				y=d[i];
				h=s*g;
				g *= c;
				z=sqrt(f*f+h*h);
				c=f/z;
				s=h/z;
				if (i1 != k1) b[i1-1]=z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				rotcol(1,n0,i1,i,v,c,s);
				d[i1]=z=sqrt(f*f+h*h);
				c=f/z;
				s=h/z;
				f=c*g+s*y;
				x=c*y-s*g;
				rotcol(1,m,i1,i,u,c,s);
			}
			b[n1]=f;
			d[n]=x;
		}
	} while (n > 0);
	em[3]=bmax;
	em[5]=count;
	em[7]=rnk;
	return n;
}
void free_real_vector(v, l)
double *v; 
int l;
{
    free((char *) (v+l));
}



void free_real_matrix(m, lr, ur, lc)
double **m;
int lr;
int ur;
int lc;
{
    int i;
    for (i=ur; i>=lr; i--) free((char *) (m[i]+lc));
    free((char *) (m+lr));
}



double *allocate_real_vector(l, u)
int l;
int u;
{
    double *p;

    p = (double *)malloc((unsigned) (u-l+1)*sizeof(double));
    if (!p) {
	fprintf(stderr, "Memory allocation failure in allocate_real_vector\n");
	exit(1);
    }
   
    return(p-l);
}



double **allocate_real_matrix(lr, ur, lc, uc)
int lr;
int ur;
int lc;
int uc;
{
    int i;
    double **p;


    p = (double **)malloc((unsigned) (ur-lr+1)*sizeof(double *));
    if (!p) {
        fprintf(stderr, "Memory allocation failure in allocate_real_matrix\n");
        exit(1);
    }
    
    p -= lr;

    for (i=lr; i<=ur; i++){
        p[i] = (double *)malloc((unsigned) (uc-lc+1)*sizeof(double));
        if (!p) {
            fprintf(stderr, "Memory allocation failure in allocate_real_matrix\n");
            exit(1);
        }
	p[i] -= lc;
    }

    return p;
}

    
void rotcol(l, u, i, j, a, c, s)
int l;
int u;
int i;
int j;
double **a;
double c;
double s;
{
	double x, y;

	for (; l<=u; l++) {
		x=a[l][i];
		y=a[l][j];
		a[l][i]=x*c+y*s;
		a[l][j]=y*c-x*s;
	}
}


void psttfmmat(a, n, v, b)
double **a;
int n;
double **v;
double b[];
{
	double matmat();
	void elmcol();
	int i,i1,j;
	double h;

	i1=n;
	v[n][n]=1.0;
	for (i=n-1; i>=1; i--) {
		h=b[i]*a[i][i1];
		if (h < 0.0) {
			for (j=i1; j<=n; j++) v[j][i]=a[i][j]/h;
			for (j=i1; j<=n; j++)
				elmcol(i1,n,j,i,v,v,matmat(i1,n,i,j,a,v));
		}
		for (j=i1; j<=n; j++) v[i][j]=v[j][i]=0.0;
		v[i][i]=1.0;
		i1=i;
	}
}
void pretfmmat(a, m, n, d)
double **a;
int m;
int n;
double d[];
{
	double tammat();
	void elmcol();
	int i,i1,j;
	double g,h;

	for (i=n; i>=1; i--) {
		i1=i+1;
		g=d[i];
		h=g*a[i][i];
		for (j=i1; j<=n; j++) a[i][j]=0.0;
		if (h < 0.0) {
			for (j=i1; j<=n; j++)
				elmcol(i,m,j,i,a,a,tammat(i1,m,i,j,a,a)/h);
			for (j=i; j<=m; j++) a[j][i] /= g;
		} else
			for (j=i; j<=m; j++) a[j][i]=0.0;
		a[i][i] += 1.0;
	}
}

void hshreabid(a, m, n, d, b, em)
double **a;
int m;
int n;
double d[];
double b[];
double em[];
{
	double tammat();
	double mattam();
	void elmcol();
	void elmrow();
	int i,j,i1;
	double norm,machtol,w,s,f,g,h;

	norm=0.0;
	for (i=1; i<=m; i++) {
		w=0.0;
		for (j=1; j<=n; j++) w += fabs(a[i][j]);
		if (w > norm) norm=w;
	}
	machtol=em[0]*norm;
	em[1]=norm;
	for (i=1; i<=n; i++) {
		i1=i+1;
		s=tammat(i1,m,i,i,a,a);
		if (s < machtol)
			d[i]=a[i][i];
		else {
			f=a[i][i];
			s += f*f;
			d[i] = g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h=f*g-s;
			a[i][i]=f-g;
			for (j=i1; j<=n; j++)
				elmcol(i,m,j,i,a,a,tammat(i,m,i,j,a,a)/h);
		}
		if (i < n) {
			s=mattam(i1+1,n,i,i,a,a);
			if (s < machtol)
				b[i]=a[i][i1];
			else {
				f=a[i][i1];
				s += f*f;
				b[i] = g = (f < 0.0) ? sqrt(s) : -sqrt(s);
				h=f*g-s;
				a[i][i1]=f-g;
				for (j=i1; j<=m; j++)
					elmrow(i1,n,j,i,a,a,mattam(i1,n,i,j,a,a)/h);
			}
		}
	}
}



void elmcol(l, u, i, j, a, b, x)
int l;
int u;
int i;
int j;
double **a;
double **b;
double x;
{
	for (; l<=u; l++) a[l][i] += b[l][j]*x;
}




void elmrow(l, u, i, j, a, b, x)
int l;
int u;
int i;
int j;
double **a;
double **b;
double x;
{
	for (; l<=u; l++) a[i][l] += b[j][l]*x;
}



double matmat(l, u, i, j, a, b)
int l;
int u;
int i;
int j;
double **a;
double **b;
{
	int k;
	double s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[i][k]*b[k][j];
	return (s);
}
