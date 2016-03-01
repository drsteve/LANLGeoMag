/*! \file Lgm_DFI_RBF.c
 *  \brief Routines to perform Divergence-Free-Interpolation of vector field data (for example, B defined on meshes).
 *
 *
 *  These routines use Divergence-Free matrix-valued Radial Basis Functions (RBF) in order to
 *  interpolate a vector field from a set of scattered data points. A
 *  vector-valued interpolant \f$\vec{s}\f$ at position vector \f$\vec{x}\f$ is
 *  computed  as follows;
 *
 *      \f[
 *          \vec{s}(\vec{x}) = \sum_{j=1}^{N} \Phi( \| \vec{x} - \vec{x}_j \| ) \vec{c}_j
 *      \f]
 *
 *  where \f$ \Phi \f$ is a matrix-valued RBF;
 *
 *      \f[
 *          \Phi = \left( \begin{array}{ccc}
 *                          \Phi_{11} & \Phi_{12} & \Phi_{13} \\
 *                          \Phi_{21} & \Phi_{22} & \Phi_{23} \\
 *                          \Phi_{31} & \Phi_{32} & \Phi_{33}
 *                        \end{array} \right).
 *      \f]
 *
 *  and \f$\vec{c}\f$ is a vector-valued weigthing factor. In order to achieve
 *  a divergence-free interpolation, the matrix-valued RBF can be computed from
 *  the following;
 *
 *      \f[
 *          \Phi = \{\nabla\nabla^T - \nabla^2 I\}\: \psi( \| \vec{x} - \vec{x}_j \| )
 *      \f]
 *
 *  where \f$ \psi( \| \vec{x} - \vec{x}_j \| )\f$ is a scalar-valued RBF
 *  (e.g., see {\em "Lowitzsch S., 2002, PhD thesis, Texas A&M University",
 *  "Lowitzsch S.  Error estimates for matrix-valued radial basis function
 *  interpolation, J. Approx. Theory, 137 (2005) 238 – 249.", "McNally, C. P.,
 *  Divergence-free interpolation of vector fields from point values - exact
 *  \f$\nabla\cdot B = 0\f$ in numerical simulations, Monthly Notices of the
 *  Royal Astronomical Society, Volume 413, Issue 1, pp. L76-L80."}).
 *
 *  In matrix form the operator \f$\nabla\nabla^T - \nabla^2 I\f$ is;
 *
 *
 *      \f{eqnarray*}
 *            {\nabla\nabla^T - \nabla^2 I} & = & { \left( \begin{array}{ccc}
 *                                                            \partial_x^2         & \partial_x\partial_y &  \partial_x\partial_z \\
 *                                                            \partial_y\partial_x & \partial_y^2         & \partial_y\partial_z  \\
 *                                                            \partial_z\partial_x & \partial_z\partial_y &          \partial_z^2  \end{array} \right)
 *                                                - \left( \begin{array}{ccc}
 *                                                             \partial_x^2 + \partial_y^2 + \partial_z^2 & 0 &0 \\
 *                                                             0 & \partial_x^2 + \partial_y^2 + \partial_z^2        & 0  \\
 *                                                             0 & 0  & \partial_x^2 + \partial_y^2 + \partial_z^2  \end{array} \right) } \\
 *                           &=& {\left( \begin{array}{ccc}
 *                                         -\partial_y^2 - \partial_z^2 & \partial_x\partial_y &  \partial_x\partial_z \\
 *                                         \partial_y\partial_x &  -\partial_x^2 - \partial_z^2        & \partial_y\partial_z  \\
 *                                         \partial_z\partial_x & \partial_z\partial_y &  -\partial_x^2 - \partial_y^2    \end{array} \right)}
 *      \f}
 *
 *
 *  If we take the \f$\psi\f$ to be a scalar guassian RBF,
 *
 *
 *      \f[
 *          \psi(r) = e^{-\epsilon r^2}
 *      \f]
 *
 *  where \f$ r^2 = x^2 + y^2 + z^2 \f$ and \f$\epsilon > 0\f$. This gives,
 *
 *      \f[
 *          \begin{array}{rclcl}
 *          \Phi_{00} & = & -\partial_y^2 - \partial_z^2 & = & \left(4\epsilon - 4\epsilon^2 \left(y^2+z^2\right)\right) \psi\\
 *          \Phi_{01} & = & \partial_x\partial_y         & = &  4 \epsilon^2 x y\: \psi\\
 *          \Phi_{02} & = & \partial_x\partial_z         & = &  4 \epsilon^2 x z\: \psi\\
 *          \Phi_{10} & = & \partial_y\partial_x         & = &  \Phi_{01}\\
 *          \Phi_{11} & = & -\partial_x^2 - \partial_z^2 & = & \left(4\epsilon - 4\epsilon^2 \left(x^2+z^2\right)\right) \psi\\
 *          \Phi_{12} & = & \partial_y\partial_z         & = &  4 \epsilon^2 y z\: \psi\\
 *          \Phi_{20} & = & \partial_z\partial_x         & = &  \Phi_{02}\\
 *          \Phi_{21} & = & \partial_z\partial_y         & = &  \Phi_{12}\\
 *          \Phi_{22} & = & -\partial_x^2 - \partial_y^2 & = & \left(4\epsilon - 4\epsilon^2 \left(x^2+y^2\right)\right) \psi
 *          \end{array}
 *      \f]
 *
 *  If we take the \f$\psi\f$ to be the multi-quadric RBF,
 *
 *
 *      \f[
 *          \psi(r) = \sqrt{1 + e r^2}
 *      \f]
 *
 *  we get;
 *
 *
 *      \f[
 *          \begin{array}{rclcl}
 *          \Phi_{00} & = & -\partial_y^2 - \partial_z^2 & = & -( \epsilon^2 (r^2+x^2) + 2\epsilon) \psi^{-3} \\
 *          \Phi_{01} & = & \partial_x\partial_y         & = & -e^2 x y \psi^{-3} \\
 *          \Phi_{02} & = & \partial_x\partial_z         & = & -e^2 x z \psi^{-3} \\
 *          \Phi_{10} & = & \partial_y\partial_x         & = &  \Phi_{01}\\
 *          \Phi_{11} & = & -\partial_x^2 - \partial_z^2 & = & -( \epsilon^2 (r^2+y^2) + 2\epsilon) \psi^{-3} \\
 *          \Phi_{12} & = & \partial_y\partial_z         & = & -e^2 y z \psi^{-3} \\
 *          \Phi_{20} & = & \partial_z\partial_x         & = &  \Phi_{02}\\
 *          \Phi_{21} & = & \partial_z\partial_y         & = &  \Phi_{12}\\
 *          \Phi_{22} & = & -\partial_x^2 - \partial_y^2 & = & -( \epsilon^2 (r^2+z^2) + 2\epsilon) \psi^{-3} 
 *          \end{array}
 *      \f]
 *
 *
 *
 *  To find the weighting vectors \f$\vec{c}_j\f$, we make sure that the
 *  interpolation agrees with all the data that we have. This leads to a linear
 *  set of equations to solve; 
 *
 *
 *      \f[
 *          {\bf A} {\bf c} = {\bf d}
 *      \f]
 *
 *  or;
 *
 *
 *      \f[ 
 *          \left( \begin{array}{cccc}
 *
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{0,0}
 *                  &
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{0,1}
 *                  &
 *                  \cdots
 *                  &
 *                  \left( \begin{array}{lll}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{0,n-1} \\
 *
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{1,0}
 *                  &
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{1,1}
 *                  &
 *                  \cdots
 *                  &
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{1,n-1} \\
 *
 *                  \vdots & \vdots & \ddots & \vdots \\
 *
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{n-1,0}
 *                  &
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{n-1,1}
 *                  &
 *                  \cdots
 *                  &
 *                  \left( \begin{array}{ccc}
 *                      \Phi_{00} & \Phi_{01} & \Phi_{02} \\
 *                      \Phi_{10} & \Phi_{11} & \Phi_{12} \\
 *                      \Phi_{20} & \Phi_{21} & \Phi_{22} \\
 *                  \end{array} \right)_{n-1,n-1} \\
 *
 *
 *           \end{array} \right)
 *          \left( \begin{array}{c}
 *                  c_{0,x} \\
 *                  c_{0,y} \\
 *                  c_{0,z} \\\noalign{\medskip}
 *                  c_{1,x} \\
 *                  c_{1,y} \\
 *                  c_{1,z} \\
 *                  \vdots  \\
 *                  c_{n-1,x} \\
 *                  c_{n-1,y} \\
 *                  c_{n-1,z} \\
 *           \end{array} \right)
 *          = 
 *          \left( \begin{array}{c}
 *                  d_{0,x} \\
 *                  d_{0,y} \\
 *                  d_{0,z} \\\noalign{\medskip}
 *                  d_{1,x} \\
 *                  d_{1,y} \\
 *                  d_{1,z} \\
 *                  \vdots  \\
 *                  d_{n-1,x} \\
 *                  d_{n-1,y} \\
 *                  d_{n-1,z} \\
 *           \end{array} \right)
 *  
 *      \f]
 *
 *  where the subscripts \f$i,j\f$ on each sub-matrix of A means ''evaluate the
 *  matrix elements at \f$(x,y,z) = (x_i-x_j, y_i-y_j, z_i-z_j)\f$''. Note that
 *  each submatrix of \f${\bf A}\f$ is symmetric, and the entire matrix should
 *  be as well because (for example);
 *
 *      \f[
 *          \Phi_{02}(x_i-x_j, y_i-y_j, z_i-z_j) = \Phi_{02}(x_j-x_i, y_j-y_i, z_j-z_i) = \Phi_{20}(x_j-x_i, y_j-y_i, z_j-z_i)
 *      \f]
 *
 *  We solve this system using Cholesky decomposition (which is faster than LU
 *  decomposition) in order to get the weigthing vectors \f$c\f$. 
 *
 *
 */

#include "Lgm/Lgm_RBF.h"

//#define LGM_DFI_RBF_SOLVER  LGM_CHOLESKY_DECOMP
#define LGM_DFI_RBF_SOLVER  LGM_PLU_DECOMP


/** Given \f$\vec{v} = (x, y, z)\f$ and \f$\vec{v}_0 = (x_0, y_0, z_0)\f$ this
 *  routine computes the matrix elements \f$\Phi_{ij}(x-x_0, y-y_0, z-z_0)\f$,
 *  where the scalar RBF is taken to be the guassian;
 *
 *  \f[
 *      \psi = \exp{ -\epsilon r^2 } = (\exp{ -\epsilon (x^2+y^2+z^2) }
 *  \f]. 
 *
 *  \param[in]          v  -   input target vector
 *  \param[in]         v0  -   input reference vector
 *  \param[in]        eps  -   smoothing factor in scalar RBF
 *  \param[in,out]    Phi  -   3x3 matrix-value RBF.
 *  \param[in]        rbf  -   pointer to structure containing info for RBF interpolation. 
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    January 24, 2012
 *
 *
 */
void    Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double Phi[3][3], Lgm_DFI_RBF_Info *rbf  ) {

    double  psi, psi2, psi3, psi5, x, y, z, x2, y2, z2, r2, f, g;
    double  eps, Eps2, TwoEps;

    eps = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        f = 4.0*eps; g = f*eps;
        psi = exp( -eps*r2 );

        Phi[0][0] = ( f - g*(y2 + z2) ) * psi;
        Phi[0][1] = g*x*y*psi;
        Phi[0][2] = g*x*z*psi;

        Phi[1][0] = Phi[0][1];
        Phi[1][1] = ( f - g*(x2 + z2) ) * psi;
        Phi[1][2] = g*y*z*psi;

        Phi[2][0] = Phi[0][2];
        Phi[2][1] = Phi[1][2];
        Phi[2][2] = ( f - g*(x2 + y2) ) * psi;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        Eps2   = eps*eps;
        TwoEps = 2.0*eps;
        f      = 1.0 + eps*r2;
        psi    = sqrt( f );
        psi2   = f; psi3 = psi*psi2;
        g      = 1.0/psi3;

        Phi[0][0] = -( Eps2*(r2+x2) + TwoEps) *g;
        Phi[0][1] = -Eps2*x*y *g;
        Phi[0][2] = -Eps2*x*z *g;


        Phi[1][0] = Phi[0][1];
        Phi[1][1] = -( Eps2*(r2+y2) + TwoEps) *g;
        Phi[1][2] = -Eps2*y*z *g;

        Phi[2][0] = Phi[0][2];
        Phi[2][1] = Phi[1][2];
        Phi[2][2] = -( Eps2*(r2+z2) + TwoEps) *g;
        
    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        Eps2   = eps*eps;
        TwoEps = 2.0*eps;
        f      = 1.0 + eps*r2;
        psi    = 1.0/sqrt( f );
        psi2   = f; psi3 = psi*psi2; psi5 = psi3*psi2;
        g      = 1.0/psi5;

        Phi[0][0] = -( Eps2*(-2.0*x2+y2+z2) - TwoEps) *g;
        Phi[0][1] = 3.0*Eps2*x*y *g;
        Phi[0][2] = 3.0*Eps2*x*z *g;


        Phi[1][0] = Phi[0][1];
        Phi[1][1] = -( Eps2*(x2-2.0*y2+z2) - TwoEps) *g;
        Phi[1][2] = 3.0*Eps2*y*z *g;

        Phi[2][0] = Phi[0][2];
        Phi[2][1] = Phi[1][2];
        Phi[2][2] = -( Eps2*(x2+y2-2.0*z2) - TwoEps) *g;
        

    } else {
        printf("Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    /*
    int i, j;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            printf("%15g ", Phi[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    */

    return;


}


void    Lgm_DFI_RBF_dPhi_dx( Lgm_Vector *v, Lgm_Vector *v0, double dPdx[3][3], Lgm_DFI_RBF_Info *rbf  ) {

    double  psi, psi2, psi3, x, y, z, x2, y2, z2, r2, f, f2, g;
    double  e, e2, e3;

    e = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi = exp( -e*r2 );

        dPdx[0][0] = 8.0*e*e*x*( e*(y2 + z2) - 1.0) * psi;
        dPdx[0][1] = 4.0*e*e*y*(1.0 - 2.0*e*x2) * psi;
        dPdx[0][2] = 4.0*e*e*z*(1.0 - 2.0*e*x2) * psi;

        dPdx[1][0] = dPdx[0][1];
        dPdx[1][1] = 8.0*e*e*x*( e*(x2 + z2) - 2.0 ) * psi;
        dPdx[1][2] = -8.0*e*e*e*x*y*z * psi;

        dPdx[2][0] = dPdx[0][2];
        dPdx[2][1] = dPdx[1][2];
        dPdx[2][2] = 8.0*e*e*x*( e*(x2 + y2) - 2.0 ) * psi;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f;
        psi = sqrt( f );
        g   = f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdx[0][0] = e2*x*(2.0 + e*(2.0*x2-y2-z2)) / g;
        dPdx[0][1] = -( e2*y*(1.0+e*(y2+z2-2.0*x2)) ) / g;
        dPdx[0][2] = -( e2*z*(1.0+e*(y2+z2-2.0*x2)) ) / g;

        dPdx[1][0] = dPdx[0][1];
        dPdx[1][1] = e2*x*(4.0 + e*(x2+4.0*y2+z2)) / g;
        dPdx[1][2] = 3.0*e3*x*y*z/ g;

        dPdx[2][0] = dPdx[0][2];
        dPdx[2][1] = dPdx[1][2];
        dPdx[2][2] = e2*x*(4.0 + e*(x2+y2+4.0*z2)) / g;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f;
        psi = sqrt( f );
        g   = f*f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdx[0][0] = 3.0*e2*x*( -2.0 + e*(-2.0*x2+3.0*(y2+z2)) ) / g;
        dPdx[0][1] = 3.0*e2*y*( 1.0 + e*(-4.0*x2+y2+z2) ) / g;
        dPdx[0][2] = 3.0*e2*z*( 1.0 + e*(-4.0*x2+y2+z2) ) / g;

        dPdx[1][0] = dPdx[0][1];
        dPdx[1][1] = 3.0*e2*x*( -4.0 + e*(x2-4.0*y2+z2) ) / g;
        dPdx[1][2] = -15.0*e3*x*y*z / g;

        dPdx[2][0] = dPdx[0][2];
        dPdx[2][1] = dPdx[1][2];
        dPdx[2][2] = 3.0*e2*x*( -4.0 + e*(x2+y2-4.0*z2) ) / g;

    } else {
        printf("Lgm_DFI_RBF_dPhi_dx: Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    return;

}


void    Lgm_DFI_RBF_dPhi_dy( Lgm_Vector *v, Lgm_Vector *v0, double dPdy[3][3], Lgm_DFI_RBF_Info *rbf  ) {

    double  psi, psi2, psi3, x, y, z, x2, y2, z2, r2, f, f2, f3, f5, g;
    double  e, e2, e3;

    e = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi = exp( -e*r2 );

        dPdy[0][0] = 8.0*e*e*y*( e*(y2 + z2) - 2.0) * psi;
        dPdy[0][1] = 4.0*e*e*x*(1.0 - 2.0*e*y2) * psi;
        dPdy[0][2] = -8.0*e*e*e*x*y*z * psi;

        dPdy[1][0] = dPdy[0][1];
        dPdy[1][1] = 8.0*e*e*y*( e*(x2 + z2) - 1.0 ) * psi;
        dPdy[1][2] = 4.0*e*e*z*(1.0 - 2.0*e*y2) * psi;

        dPdy[2][0] = dPdy[0][2];
        dPdy[2][1] = dPdy[1][2];
        dPdy[2][2] = 8.0*e*e*y*( e*(x2 + y2) - 2.0 ) * psi;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f;
        psi = sqrt( f );
        g   = f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdy[0][0] = e2*y*(4.0 + e*(4.0*x2+y2+z2)) / g;
        dPdy[0][1] = -( e2*x*(1.0+e*(x2-2.0*y2+z2)) ) / g;
        dPdy[0][2] = 3.0*e3*x*y*z / g;

        dPdy[1][0] = dPdy[0][1];
        dPdy[1][1] = e2*y*(2.0 - e*(x2-2.0*y2+z2)) / g;
        dPdy[1][2] = -( e2*z*(1.0+e*(x2-2.0*y2+z2)) ) / g;

        dPdy[2][0] = dPdy[0][2];
        dPdy[2][1] = dPdy[1][2];
        dPdy[2][2] = e2*y*(4.0 + e*(x2+y2+4.0*z2)) / g;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f;
        psi = sqrt( f );
        g   = f*f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdy[0][0] = 3.0*e2*y*( -4.0 + e*(-4.0*x2+y2+z2) ) / g;
        dPdy[0][1] = 3.0*e2*x*( 1.0 + e*(x2-4.0*y2+z2) ) / g;
        dPdy[0][2] = -15.0*e3*x*y*z / g;

        dPdy[1][0] = dPdy[0][1];
        dPdy[1][1] = 3.0*e2*y*( -2.0 + e*(3.0*(x2+z2)-2.0*y2) ) / g;
        dPdy[1][2] = 3.0*e2*z*( 1.0 + e*(x2-4.0*y2+z2) ) / g;

        dPdy[2][0] = dPdy[0][2];
        dPdy[2][1] = dPdy[1][2];
        dPdy[2][2] = 3.0*e2*y*( -4.0 + e*(x2+y2-4.0*z2) ) / g;

    } else {
        printf("Lgm_DFI_RBF_dPhi_dy: Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    return;

}


void    Lgm_DFI_RBF_dPhi_dz( Lgm_Vector *v, Lgm_Vector *v0, double dPdz[3][3], Lgm_DFI_RBF_Info *rbf  ) {

    double  psi, psi2, psi3, x, y, z, x2, y2, z2, r2, f, f2, f3, f5, g;
    double  e, e2, e3;

    e = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi = exp( -e*r2 );

        dPdz[0][0] = 8.0*e*e*z*( e*(y2 + z2) - 2.0) * psi;
        dPdz[0][1] = -8.0*e*e*e*x*y*z * psi;
        dPdz[0][2] = 4.0*e*e*x*(1.0 - 2.0*e*z2) * psi;

        dPdz[1][0] = dPdz[0][1];
        dPdz[1][1] = 8.0*e*e*z*( e*(x2 + z2) - 2.0 ) * psi;
        dPdz[1][2] = 4.0*e*e*y*(1.0 - 2.0*e*z2) * psi;

        dPdz[2][0] = dPdz[0][2];
        dPdz[2][1] = dPdz[1][2];
        dPdz[2][2] = 8.0*e*e*z*( e*(x2 + y2) - 1.0 ) * psi;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f; 
        psi = sqrt( f );
        g   = f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdz[0][0] = e2*z*(4.0 + e*(4.0*x2+y2+z2)) / g;
        dPdz[0][1] = 3.0*e3*x*y*z / g;
        dPdz[0][2] = -( e2*x*(1.0+e*(x2+y2-2.0*z2)) ) / g;

        dPdz[1][0] = dPdz[0][1];
        dPdz[1][1] = e2*z*(4.0 + e*(x2+4.0*y2+z2)) / g;
        dPdz[1][2] = -( e2*y*(1.0+e*(x2+y2-2.0*z2)) ) / g;

        dPdz[2][0] = dPdz[0][2];
        dPdz[2][1] = dPdz[1][2];
        dPdz[2][2] = e2*z*(2.0 - e*(x2+y2-2.0*z2)) / g;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2; f2 = f*f;
        psi = sqrt( f );
        g   = f*f2*psi;
        e2  = e*e; e3 = e2*e;

        dPdz[0][0] = 3.0*e2*z*( -4.0 + e*(-4.0*x2+y2+z2) ) / g;
        dPdz[0][1] = -15.0*e3*x*y*z / g;
        dPdz[0][2] = 3.0*e2*x*( 1.0 + e*(x2+y2-4.0*z2) ) / g;

        dPdz[1][0] = dPdz[0][1];
        dPdz[1][1] = 3.0*e2*z*( -4.0 + e*(x2-4.0*y2+z2) ) / g;
        dPdz[1][2] = 3.0*e2*y*( 1.0 + e*(x2+y2-4.0*z2) ) / g;

        dPdz[2][0] = dPdz[0][2];
        dPdz[2][1] = dPdz[1][2];
        dPdz[2][2] = 3.0*e2*z*( -2.0 + e*(3.0*(x2+y2)-2.0*z2) ) / g;

    } else {
        printf("Lgm_DFI_RBF_dPhi_dz: Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    return;

}




/** From a vector-field dataset, compute the vector-valued weighting factors,
 *  \f$\vec{c}_j\f$. Info is returned in the rbf structure.
 *
 *
 *  \param[in]                       v   -   pointer to an array of position vectors.
 *  \param[in]                       B   -   pointer to array of corresponding field vectors.
 *  \param[in]                       n   -   number of (v, B) pairs defined.
 *  \param[in]                     eps   -   smoothing factor in scalar RBF.
 *  \param[in]      RadialBasisFunction  -   RBF to use. Can be LGM_RBF_GAUSSIAN, LGM_RBF_MULTIQUADRIC
 *
 *  \return  pointer to structure containing info for RBF interpolation. User
 *           is responsible for freeing with Lgm_DFI_RBF_Free().
 *
 *  \author  M. G. Henderson
 *  date    January 24, 2012
 *
 *
 */
Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( unsigned long int *I_data, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction ) {

    int              i, j, ii, jj, p, q, n3, s;
    double           *d, **a, Phi[3][3], val;
    gsl_matrix       *A, *V;
    gsl_vector       *D, *c, *S, *Work;
    Lgm_DFI_RBF_Info *rbf;

    n3 = 3*n;
    A = gsl_matrix_calloc( n3, n3 );
    c = gsl_vector_alloc( n3 );
    D = gsl_vector_calloc( n3 );


    /*
     * Save info needed to do an evaluation.
     */
    rbf = ( Lgm_DFI_RBF_Info *)calloc( 1, sizeof(*rbf) );
    rbf->RadialBasisFunction = RadialBasisFunction;
    rbf->eps = eps;
    rbf->n   = n;
    rbf->n3  = n3;
    LGM_ARRAY_1D( rbf->LookUpKey, n, unsigned long int);
    LGM_ARRAY_1D( rbf->v, n, Lgm_Vector);
    LGM_ARRAY_1D( rbf->c, n, Lgm_Vector);
    for ( i=0; i<n; i++ ) {
        rbf->LookUpKey[i] = I_data[i];
        rbf->v[i] = v[i];
    }
    // This subtraction doesntm seem to work out very well...
//    rbf->Bx0 = B[0].x;
//    rbf->By0 = B[0].y;
//    rbf->Bz0 = B[0].z;

double Bbkg;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].x; rbf->Bx0 = Bbkg/(double)n;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].y; rbf->By0 = Bbkg/(double)n;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].z; rbf->Bz0 = Bbkg/(double)n;
    rbf->Bx0 = 0.0;
    rbf->By0 = 0.0;
    rbf->Bz0 = 0.0;
    
    /*
     * Fill d array. (Subtract off the field at the nearest point v[0] -- See
     * McNally [2011].) We add this field back on later.
     */
    for (i=0; i<n; i++){
        gsl_vector_set( D, 3*i+0, B[i].x - rbf->Bx0 );
        gsl_vector_set( D, 3*i+1, B[i].y - rbf->By0 );
        gsl_vector_set( D, 3*i+2, B[i].z - rbf->Bz0 );
    }


    /*
     *                                             [  row0  ]
     * Fill A matrix. In C, order is A[row][col] = [  row1  ]
     *                                             [  row2  ]
     */
    for ( i=0; i<n; i++ ) { // locate start row for subarray
        ii = 3*i;

        for ( j=0; j<n; j++ ) { // locate start column for subarray
            jj = 3*j;

            // Get Phi( v_i - v_j )
            Lgm_DFI_RBF_Phi( &v[i], &v[j], Phi, rbf );

            for ( p=0; p<3; p++ ){ // subarray row
                for ( q=0; q<3; q++ ){  // subarray column
                    gsl_matrix_set( A, ii+p, jj+q, Phi[p][q] );
                }
            }


        }

    }

    /*
    for (i=0; i<n; i++ ) {
        printf("v%02d = %8g %8g %8g   B%02d = %8g %8g %8g\n", i, v[i].x, v[i].y, v[i].z, i, B[i].x, B[i].y, B[i].z );
    }
    for (i=0; i<n3; i++){
        for (j=0; j<n3; j++){
            printf("%8g ", gsl_matrix_get(A, i, j ) );
        }
        printf("\n");
    }
    */




    /*
     * Now we need to solve the system of equation;
     *
     *      d = ac
     *
     *  for c.
     *
     *  First create gsl_vector and gsl_matrix views of the d and A arrays.
     *  Then compute Cholesky decomposition of the a array. Then solve the
     *  system to get c.
     *
     */
    if ( LGM_DFI_RBF_SOLVER == LGM_CHOLESKY_DECOMP ){
        gsl_linalg_cholesky_decomp( A );
        gsl_linalg_cholesky_solve( A, D, c );
    } else if ( LGM_DFI_RBF_SOLVER == LGM_PLU_DECOMP ){
        gsl_permutation *P = gsl_permutation_alloc( n3 );
        gsl_linalg_LU_decomp( A, P, &s );
        gsl_linalg_LU_solve( A, P, D, c );
        gsl_permutation_free( P );
    } else if ( LGM_DFI_RBF_SOLVER == LGM_SVD ){
        V    = gsl_matrix_calloc( n3, n3 );
        S    = gsl_vector_alloc( n3 );
        Work = gsl_vector_alloc( n3 );
        gsl_linalg_SV_decomp( A, V, S, Work );
        gsl_linalg_SV_solve( A, V, S, D, c );
        gsl_vector_free( Work );
        gsl_vector_free( S );
        gsl_matrix_free( V );
    }

    for (i=0; i<n; i++){
        rbf->c[i].x = gsl_vector_get( c, 3*i+0 );
        rbf->c[i].y = gsl_vector_get( c, 3*i+1 );
        rbf->c[i].z = gsl_vector_get( c, 3*i+2 );
    }


    
    gsl_vector_free( D );
    gsl_vector_free( c );
    gsl_matrix_free( A );

    return( rbf );

}

/** Free a previously allocated Lgm_DFI_RBF_Info structure.
 *
 *
 *  \param[in]     rbf   -   pointer to structure containing info for RBF interpolation. 
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    January 24, 2012
 *
 *
 */
void    Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf ) {
    LGM_ARRAY_1D_FREE( rbf->LookUpKey );
    LGM_ARRAY_1D_FREE( rbf->v );
    LGM_ARRAY_1D_FREE( rbf->c );
    free( rbf );
    return;
}







/**  Compute Divergence Free interpolant at the specified position vector. The
 *   weights given in \f$\vec{c}\f$ must have been pre-computed with
 *   Lgm_DFI_RBF_Init().
 *
 *
 *  \param[in]        v   -   position vector to compute B at.
 *  \param[out]       B   -   interpolated value of B at v.
 *  \param[out]     rbf   -   pointer to initialized Lgm_DFI_RBF_Info structure.
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    January 24, 2012
 *
 *
 */
void    Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf ) {

    int         j;
    Lgm_Vector  W;
    double      Phi[3][3];


    B->x = B->y = B->z = 0.0;
    for ( j=0; j<rbf->n; j++ ){

        Lgm_DFI_RBF_Phi( v, &rbf->v[j], Phi, rbf );
        Lgm_MatTimesVec( Phi, &rbf->c[j], &W );
        //printf("rbf->n = %d   c = %g %g %g   W = %g %g %g\n", rbf->n, rbf->c[j].x, rbf->c[j].y, rbf->c[j].z, W.x, W.y, W.z);
        B->x += W.x; B->y += W.y; B->z += W.z;

    }

    // Add the subtracted "background" back in (that we subtracted prior to the interp).
    B->x += rbf->Bx0;
    B->y += rbf->By0;
    B->z += rbf->Bz0;
    //printf("rbf->Bx0, rbf->By0, rbf->Bz0 = %g %g %g\n", rbf->Bx0, rbf->By0, rbf->Bz0 );
    //printf("v = %g %g %g   B = %g %g %g\n", v->x, v->y, v->z, B->x, B->y, B->z );


    return;
    
}





/**  Compute Divergence Free interpolant at the specified position vector. The
 *   weights given in \f$\vec{c}\f$ must have been pre-computed with
 *   Lgm_DFI_RBF_Init().
 *
 *
 *  \param[in]         v   -   position vector to compute B at.
 *  \param[out]     dBdx   -   interpolated value of dB/dx at v.
 *  \param[out]     dBdy   -   interpolated value of dB/dy at v.
 *  \param[out]     dBdz   -   interpolated value of dB/dz at v.
 *  \param[out]      rbf   -   pointer to initialized Lgm_DFI_RBF_Info structure.
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    July 1, 2015
 *
 *
 */
void    Lgm_DFI_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_DFI_RBF_Info *rbf ) {

    int         j;
    Lgm_Vector  W;
    double      P[3][3];


    dBdx->x = dBdx->y = dBdx->z = 0.0;
    dBdy->x = dBdy->y = dBdy->z = 0.0;
    dBdz->x = dBdz->y = dBdz->z = 0.0;
    for ( j=0; j<rbf->n; j++ ){

        Lgm_DFI_RBF_dPhi_dx( v, &rbf->v[j], P, rbf );
        Lgm_MatTimesVec( P, &rbf->c[j], &W );
        dBdx->x += W.x; dBdx->y += W.y; dBdx->z += W.z;

        Lgm_DFI_RBF_dPhi_dy( v, &rbf->v[j], P, rbf );
        Lgm_MatTimesVec( P, &rbf->c[j], &W );
        dBdy->x += W.x; dBdy->y += W.y; dBdy->z += W.z;

        Lgm_DFI_RBF_dPhi_dz( v, &rbf->v[j], P, rbf );
        Lgm_MatTimesVec( P, &rbf->c[j], &W );
        dBdz->x += W.x; dBdz->y += W.y; dBdz->z += W.z;

    }


    return;
    
}


