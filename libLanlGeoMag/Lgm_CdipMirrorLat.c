/*! \file Lgm_CdipMirrorLat.c
 *
 *  \brief Routines for computing the centered dipole lattitude of the mirror point in a centered dipole field model.
 *
 */
#include "Lgm/Lgm_MagModelInfo.h"



/**
 *  \brief
 *      Function used to find the mirror latitude in a pure centered dipole field.
 *
 *  \details
 *      Computes,
 *
 *          \f[ F(x) = x^6 + 3s^4x - 4s^4 \f]
 *
 *      where, \f$s^4\f$ is equal to \f$\sin^4\alpha_\circ\f$
 *      (\f$\alpha_\circ\f$ is the equatorial pitch angle).  The value of
 *      \f$x\f$ will be equal to \f$\cos^2\lambda_m\f$ when \f$F(x)=0\f$. 
 *
 *
 *      \param[in]      x       A candidate value of \f$\cos^2\lambda_m\f$.
 *      \param[in]      s4      The value of \f$\sin^4\alpha_\circ\f$.
 *      \param[in]      Info    Additional information that can be passed by Lgm_zBrent(), but not needed here.
 *
 *      \return         \f$ F(x) = x^6 + 3s^4x - 4s^4 \f$
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
double MirrorLatFunc( double x, double s4, void *Info ) {
    double  x2, x3;
    x2 = x*x; x3 = x2*x;
    return( x3*x3 + s4*(3.0*x - 4.0) );
}

/**
 *  \brief
 *      Computes the mirror latitude in a pure centered dipole field for a given equatorial pitch angle.
 *
 *  \details
 *      The Magnetic field strength in a dipole is given by;
 *  
 *          \f[ B = {M\over r^3_\circ} {[4-3\cos^2\lambda ]^{1/2}\over \cos^6\lambda } \f]
 *
 *      Since,
 *
 *          \f[ {\sin^2\alpha_1\over B_1 } = { \sin^2\alpha_2\over B_2 }, \f]
 *      
 *      the mirror latitude (where local pitch angle, \f$\alpha_m = 90^\circ\f$
 *      and \f$B=B_m\f$) is related to the equatorial pitch angle
 *      (\f$\alpha_\circ\f$) by;
 *      
 *          \f[ \sin^2\alpha_\circ = {\cos^6\lambda_m\over[ 4 - 3\cos^2\lambda_m ]^{1/2} }, \f]
 *
 *      Let \f$x=\cos^2\lambda_m\f$ and \f$s=\sin\alpha_\circ\f$. Then,
 *
 *          \f[ x^6 + 3s^4x - 4s^4 = 0 \f]
 *
 *      This can be solved efficiently for x using Lgm_zBrent().
 *
 *      \param[in]      SinAlpha0   Sine of equatorial pitch angle, \f$\sin(\alpha_\circ)\f$.
 *
 *      \return         Cosine of the mirror latitude, \f$\cos(\lambda_m)\f$
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
double Lgm_CdipMirrorLat( double SinAlpha0 ){

    BrentFuncInfo   f;
    double          s2, s4, x, g;

    s2 = SinAlpha0*SinAlpha0; s4 = s2*s2;
    if (s4>1.0) s4 = 1.0;
    else if ( s4<-1.0) s4 = -1.0;

    f.func = MirrorLatFunc;
    f.Val  = s4;

    Lgm_zBrent( 0.0, 1.0, -4.0*s4, 1.0-s4, &f, 1e-10, &x, &g );
    return( sqrt(x) );

}
