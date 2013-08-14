/* 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "Lgm/Lgm_MagModelInfo.h"



int Lgm_TraceToYZPlane( Lgm_Vector *u, Lgm_Vector *v, double Xtarget, double sgn0, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s, sgn, fsgn0;
    double	    Sa, Sb, B, f, r2, z2, r3, L, R, hhh;
    double      Stotal;
    Lgm_Vector	Btmp;
    Lgm_Vector	Pa, Pb, P;
    int		    done, reset=TRUE, bracketed=FALSE;
    long int    Ntotal;

//printf("\n\n\n***************************************\n");
    Pb.x = Pb.y = Pb.z = 0.0;


    Hmax = 20.0;
    Hmin = 0.01;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;


    /* 
     *  Set the start point, Pa 
     */
    Pa   = *u;
    f = Xtarget - Pa.x;
    fsgn0 = ( f<0.0) ? -1.0 : 1.0;
    Sa   = 0.0;
    Stotal = 0.0;
    Ntotal = 0;





    /*
     *  Choose a good step size to try.
     */
    Lgm_Convert_Coords( &Pa, &P, GSM_TO_SM, Info->c );
    R  = Lgm_Magnitude( &Pa );
    r2 = R*R;
    z2 = P.z*P.z;
    L  = 9e99;
    if ( fabs( r2 - z2 ) < 1e-2 ) {

        /*
         *  High  L
         */
        Htry = 2.0;

    } else {

        r3   = r2*R;
        L    = r3 / (r2 - z2);
        Htry = ( L < 4.5 ) ? 1.165*L : 4.5;

    }

    if ( Htry < 1e-4 ) Htry = 0.001;



//Htry = 0.01;


    /* 
     *   Now, bracket minimum. And do a bisection search for the minimum.
     */
    done = FALSE;
    sgn  = sgn0;
    P    = Pa;
    Sb   = -9e99;
//printf("Htry =%g   sgn = %g\n", Htry, sgn );
//printf("Xtarget = %g f = %g P = %g %g %g\n", Xtarget, f, P.x, P.y, P.z);
    while ( !done ) {

        //if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Stotal += Hdid;
        ++Ntotal;
        R  = Lgm_Magnitude( &P );
        r2 = R*R;
        r3   = r2*R;
        L    = r3 / (r2 - z2);

        f = Xtarget - P.x;
//printf("Xtarget = %g f = %g P = %g %g %g\n", Xtarget, f, P.x, P.y, P.z);

         if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) || ( s > 1000.0 ) ) {
            /*
             *  Open FL!
             */
            //v->x = v->y = v->z = 0.0;
            *v = P;
//printf("OPEN!\n\n\n");
            return(0);

	    } else if ( r2 < 1.0 ) {

	        return( -1 );

        } else if ( (Stotal > 300.0) || ( Ntotal > 1500 ) ) {

            printf("Lgm_TraceToYZPlane: %s:%d Field line too long or too many iterations: Stotal / Ntotal:  %g / %ld\n", __FILE__, __LINE__, Stotal, Ntotal );

	        return( 0 );

	    } else if ( fabs( f ) < tol ){

	        done = TRUE;

        } else if ( fsgn0*f > 0.0 ){

            Pa    = P;
            Sa   += Hdid;
            sgn   = sgn0;
            if ( bracketed ) {
                Htry = fabs( Sb - Sa )/2.0;
            } else if ( L < 4.5 ) {
                hhh = 1.165*L+1.0 - R + 1.0;
                //if ( hhh > 1e-4 ) Htry = hhh;
                if ( hhh > R-1.0 ) {
                    hhh = (R-1.0)/2.0;
                }
                

	        }
	    /*
printf("A: Sa, Sb, P = %g %g    (%g, %g, %g)   Htry = %f\n", Sa, Sb, P.x, P.y, P.z, Htry);
*/

	    } else {

            Pb  = P;
            if ( bracketed ) Sb += sgn*sgn0*Hdid;
            else	     Sb  = Sa + Hdid;
            
            bracketed = TRUE;
            sgn  = -sgn0;
            Htry = fabs( Sb - Sa )/2.0;
	    /*
printf("B: Sa, Sb, P = %g %g    (%g, %g, %g)   Htry = %f\n", Sa, Sb, P.x, P.y, P.z, Htry);
*/

	    }

    }




    /*
     *  
     */
    *v = Pb;
//printf("B: Sa, Sb, P = %g %g    (%g, %g, %g)   Htry = %f\n", Sa, Sb, P.x, P.y, P.z, Htry);


//    if ( fabs( Xtarget - Pb.x ) >  2.0*tol ) return( -1 );

    return( 1 );

}

