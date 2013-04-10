#include <stdio.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <Lgm/Lgm_DynamicMemory.h>



#define MaxCurvePoints 5000

void MakeTube(double *X, double *Y, double *Z, int NumCurvePoints, int NumCirclePoints, double TubeRadius ){

	double CurvePoint[MaxCurvePoints][3], TangentVector[MaxCurvePoints][3]; 
	double a, b, c, L, Linv, min, AngInc, Angle;
	double x[MaxCurvePoints][3], y[MaxCurvePoints][3], z[MaxCurvePoints][3], tmp[3];
	double CirclePoints[100][3], CirclePoints_p[100][3];
	double ***Mesh, ***Normals;
	double v1[3], v2[3], v3[3], n1[3], n2[3], n3[3];
	
	int     i, j, k, n;


    LGM_ARRAY_3D( Mesh,    MaxCurvePoints, 100, 3, double );
    LGM_ARRAY_3D( Normals, MaxCurvePoints, 100, 3, double );

	/*
	 *  Read in the points
	 */
	for (n=0; n<NumCurvePoints; ++n){
	    CurvePoint[n][0] = X[n];
	    CurvePoint[n][1] = Y[n];
	    CurvePoint[n][2] = Z[n];
	}


    //TubeRadius = 0.075/2.0;
    //NumCirclePoints = 24;


    //printf("NumCurvePoints = %d\n", NumCurvePoints);
    //printf("NumCirclePoints = %d\n", NumCirclePoints);
    //printf("TubeRadius = %f\n", TubeRadius);



	/*
	 *  Compute the tangent vector at each point on the curve.
 	 */

	/* for i = 0 */
	TangentVector[0][0] = 0.5*( -3.0*CurvePoint[0][0] + 4.0*CurvePoint[1][0] - CurvePoint[2][0] );
	TangentVector[0][1] = 0.5*( -3.0*CurvePoint[0][1] + 4.0*CurvePoint[1][1] - CurvePoint[2][1] );
	TangentVector[0][2] = 0.5*( -3.0*CurvePoint[0][2] + 4.0*CurvePoint[1][2] - CurvePoint[2][2] );

	/* for i = 1 */
	TangentVector[1][0] = 0.5*( CurvePoint[2][0] - CurvePoint[0][0] );
	TangentVector[1][1] = 0.5*( CurvePoint[2][1] - CurvePoint[0][1] );
	TangentVector[1][2] = 0.5*( CurvePoint[2][2] - CurvePoint[0][2] );

   
	for (i=2; i < NumCurvePoints-2; ++i){
		TangentVector[i][0] = (CurvePoint[i-2][0] - 8.0*CurvePoint[i-1][0]
			+ 8.0*CurvePoint[i+1][0] - CurvePoint[i+2][0])/12.0;
		TangentVector[i][1] = (CurvePoint[i-2][1] - 8.0*CurvePoint[i-1][1]
			+ 8.0*CurvePoint[i+1][1] - CurvePoint[i+2][1])/12.0;
		TangentVector[i][2] = (CurvePoint[i-2][2] - 8.0*CurvePoint[i-1][2]
			+ 8.0*CurvePoint[i+1][2] - CurvePoint[i+2][2])/12.0;
	}

	/* for i = NumCurvePoints-1 */
	i = NumCurvePoints-1;
	TangentVector[i][0] = 0.5*(3.0*CurvePoint[i][0] - 4.0*CurvePoint[i-1][0] + CurvePoint[i-2][0] );
	TangentVector[i][1] = 0.5*(3.0*CurvePoint[i][1] - 4.0*CurvePoint[i-1][1] + CurvePoint[i-2][1] );
	TangentVector[i][2] = 0.5*(3.0*CurvePoint[i][2] - 4.0*CurvePoint[i-1][2] + CurvePoint[i-2][2] );

	/* for i = NumCurvePoints-2 */
	i = NumCurvePoints-2;
	TangentVector[i][0] = 0.5*( CurvePoint[i+1][0] - CurvePoint[i-1][0] );
	TangentVector[i][1] = 0.5*( CurvePoint[i+1][1] - CurvePoint[i-1][1] );
	TangentVector[i][2] = 0.5*( CurvePoint[i+1][2] - CurvePoint[i-1][2] );











	/*
	 *  Normalize the tangent vectors and rename as x.
	 */
	for (i=0; i < NumCurvePoints; ++i){
		a = TangentVector[i][0], b = TangentVector[i][1], c = TangentVector[i][2];
		L = sqrt(a*a + b*b +c*c);
		x[i][0] = TangentVector[i][0]/L;
		x[i][1] = TangentVector[i][1]/L;
		x[i][2] = TangentVector[i][2]/L;
	}








	
	/*
	 *  Now construct a right-handed coordinate system at each point from the tangent 
	 *  vector there.  Its going to be arbitrarily oriented.
	 *
	 *  First, get another vector that is not parallel to the tangent vector;
	 */
	tmp[0] = 0.0, tmp[1] = 0.0, tmp[2] = 0.0;
	min = 999.0;
	k = 0;
	for (j=0; j<3; ++j){
		if (fabs(x[0][j]) < min){
			min = fabs(x[0][j]);
			//sgn = (x[0][j] < 0.0) ? -1.0 : 1.0;
			k = j;
		}
	}
	tmp[k] = 1.0;


	/*
	 *  Now take the cross product of the tangent and tmp vector to get
	 *  another vector that is perpendiculat to both.
	 */
	y[0][0] = x[0][1]*tmp[2] - x[0][2]*tmp[1];
	y[0][1] = x[0][2]*tmp[0] - x[0][0]*tmp[2];
	y[0][2] = x[0][0]*tmp[1] - x[0][1]*tmp[0];

	/*
	 *  Now complete the right handed system. z = x cross y
	 */
	z[0][0] = x[0][1]*y[0][2] - x[0][2]*y[0][1];
	z[0][1] = x[0][2]*y[0][0] - x[0][0]*y[0][2];
	z[0][2] = x[0][0]*y[0][1] - x[0][1]*y[0][0];

	/*
	 *  Now do for all the other points
	 */
	for (i=1; i < NumCurvePoints; ++i){

		/*
		 *   The problem is that we want the polygon points to be
		 *   reasonably well lined up, so that there is no radical twist
		 *   in the tube from point to point. Soooo...
		 *   Use the y vector of the previous coord sys. crossed with x of the new
		 *   as the new tmp.
		 *   Note that we cannot in general use the same tmp vector for all
		 *   points that we did for the first one, because the curve may take
		 *   turns into directions that make the tangent the same as that first
		 *   temp vector..... (and then the cross product is undefined!)
		 */
		tmp[0] = y[i-1][1]*x[i][2] - y[i-1][2]*x[i][1];
		tmp[1] = y[i-1][2]*x[i][0] - y[i-1][0]*x[i][2];
		tmp[2] = y[i-1][0]*x[i][1] - y[i-1][1]*x[i][0];

		
		y[i][0] = x[i][1]*tmp[2] - x[i][2]*tmp[1];
		y[i][1] = x[i][2]*tmp[0] - x[i][0]*tmp[2];
		y[i][2] = x[i][0]*tmp[1] - x[i][1]*tmp[0];

		z[i][0] = x[i][1]*y[i][2] - x[i][2]*y[i][1];
		z[i][1] = x[i][2]*y[i][0] - x[i][0]*y[i][2];
		z[i][2] = x[i][0]*y[i][1] - x[i][1]*y[i][0];
		
	}


	

	/*
	 *  O.K., now we have a RH coord system at each point along the curve such that
	 *  the x axis is perpendicular to the circle we want to represent. I.e. the circle
	 *  should be drawn in the yz plane in this system.
	 */
	AngInc = 2.0*M_PI/(double)NumCirclePoints;
	Angle = 0.0;
	for (j=0; j<NumCirclePoints; ++j){
		CirclePoints_p[j][0] = 0.0;
		CirclePoints_p[j][1] = TubeRadius * cos(Angle);
		CirclePoints_p[j][2] = TubeRadius * sin(Angle);
		Angle += AngInc;
	}
	CirclePoints_p[NumCirclePoints][0] = CirclePoints_p[0][0];
	CirclePoints_p[NumCirclePoints][1] = CirclePoints_p[0][1];
	CirclePoints_p[NumCirclePoints][2] = CirclePoints_p[0][2];
 




	/*
	 * Now compute the CirclePoints_p vectors relative to the original coordinate
	 * system. Dump results in the Mesh array.
	 *      Mesh[curve point number][circle pointm number][component]
	 */
	for (i=0; i < NumCurvePoints; ++i){
		for (j=0; j<NumCirclePoints+1; ++j){
			CirclePoints[j][0] = CirclePoints_p[j][0] * x[i][0]
					   + CirclePoints_p[j][1] * y[i][0]
					   + CirclePoints_p[j][2] * z[i][0];
			CirclePoints[j][1] = CirclePoints_p[j][0] * x[i][1]
					   + CirclePoints_p[j][1] * y[i][1]
					   + CirclePoints_p[j][2] * z[i][1];
			CirclePoints[j][2] = CirclePoints_p[j][0] * x[i][2]
					   + CirclePoints_p[j][1] * y[i][2]
					   + CirclePoints_p[j][2] * z[i][2];
			
			Mesh[i][j][0] = CurvePoint[i][0] + CirclePoints[j][0];
			Mesh[i][j][1] = CurvePoint[i][1] + CirclePoints[j][1];
			Mesh[i][j][2] = CurvePoint[i][2] + CirclePoints[j][2];
		}

	}
	

	/*
	 *  Now, Compute the normalized normals at each point on the Mesh.
	 */
	for (i=0; i < NumCurvePoints; ++i){
		for (j=0; j<NumCirclePoints+1; ++j){
			Normals[i][j][0] = Mesh[i][j][0] - CurvePoint[i][0];
			Normals[i][j][1] = Mesh[i][j][1] - CurvePoint[i][1];
			Normals[i][j][2] = Mesh[i][j][2] - CurvePoint[i][2];
			a = Normals[i][j][0], b = Normals[i][j][1], c = Normals[i][j][2];
			Linv = sqrt(a*a + b*b +c*c);
			Normals[i][j][0] *= Linv;
			Normals[i][j][1] *= Linv;
			Normals[i][j][2] *= Linv;
		}
	}




	 

	/*
	 *  Now, Mesh should contain all the points on the tubluar mesh we plan
	 *  to render with trinagles. And Normals should contain all the Normal vectors
	 *  to the surface. So now all we need to do is to output the mesh in a form 
	 *  that POVray (for example) can use. Use the smooth_triangle object....
	 */

//	fp = fopen("Line.pov", "w");
//	fprintf(fp, "mesh{\n");
    glBegin( GL_TRIANGLES );
	for (i=0; i < NumCurvePoints-1; ++i){
		for (j=0; j<NumCirclePoints; ++j){

			/* 1st trinagle */
            v1[0] = Mesh[i][j][0],          v1[1] = Mesh[i][j][1],          v1[2] = Mesh[i][j][2];
            n1[0] = Normals[i][j][0],       n1[1] = Normals[i][j][1],       n1[2] = Normals[i][j][2];

            v2[0] = Mesh[i][j+1][0],        v2[1] = Mesh[i][j+1][1],        v2[2] = Mesh[i][j+1][2];
            n2[0] = Normals[i][j+1][0],     n2[1] = Normals[i][j+1][1],     n2[2] = Normals[i][j+1][2];

            v3[0] = Mesh[i+1][j+1][0],      v3[1] = Mesh[i+1][j+1][1],      v3[2] = Mesh[i+1][j+1][2];
            n3[0] = Normals[i+1][j+1][0],   n3[1] = Normals[i+1][j+1][1],   n3[2] = Normals[i+1][j+1][2];



            // In openGL, the normals come first. One normal for each vertex gives
            // smoothly-shaded triangles.
            glNormal3f( n1[0], n1[1], n1[2] ); glVertex3f( v1[0], v1[1], v1[2] );  
            glNormal3f( n2[0], n2[1], n2[2] ); glVertex3f( v2[0], v2[1], v2[2] );  
            glNormal3f( n3[0], n3[1], n3[2] ); glVertex3f( v3[0], v3[1], v3[2] );  

//			fprintf(fp, "smooth_triangle { <%f, %f, %f>, <%f, %f, %f>,", v1[0], v1[1], v1[2], n1[0], n1[1], n1[2]);
//			fprintf(fp, " <%f, %f, %f>, <%f, %f, %f>,", v2[0], v2[1], v2[2], n2[0], n2[1], n2[2]);
//			fprintf(fp, " <%f, %f, %f>, <%f, %f, %f>}\n", v3[0], v3[1], v3[2], n3[0], n3[1], n3[2]);

			/* 2nd trinagle */
			v1[0] = Mesh[i][j][0],          v1[1] = Mesh[i][j][1],          v1[2] = Mesh[i][j][2];
			n1[0] = Normals[i][j][0],       n1[1] = Normals[i][j][1],       n1[2] = Normals[i][j][2];

			v2[0] = Mesh[i+1][j+1][0],      v2[1] = Mesh[i+1][j+1][1],      v2[2] = Mesh[i+1][j+1][2];
			n2[0] = Normals[i+1][j+1][0],   n2[1] = Normals[i+1][j+1][1],   n2[2] = Normals[i+1][j+1][2];

			v3[0] = Mesh[i+1][j][0],        v3[1] = Mesh[i+1][j][1],        v3[2] = Mesh[i+1][j][2];
			n3[0] = Normals[i+1][j][0],     n3[1] = Normals[i+1][j][1],     n3[2] = Normals[i+1][j][2];

            glNormal3f( n1[0], n1[1], n1[2] ); glVertex3f( v1[0], v1[1], v1[2] );  
            glNormal3f( n2[0], n2[1], n2[2] ); glVertex3f( v2[0], v2[1], v2[2] );  
            glNormal3f( n3[0], n3[1], n3[2] ); glVertex3f( v3[0], v3[1], v3[2] );  


//			fprintf(fp, "smooth_triangle { <%f, %f, %f>, <%f, %f, %f>,", v1[0], v1[1], v1[2], n1[0], n1[1], n1[2]);
//			fprintf(fp, " <%f, %f, %f>, <%f, %f, %f>,", v2[0], v2[1], v2[2], n2[0], n2[1], n2[2]);
//			fprintf(fp, " <%f, %f, %f>, <%f, %f, %f>}\n", v3[0], v3[1], v3[2], n3[0], n3[1], n3[2]);
		}
	}
    glEnd( );
//	fprintf(fp, "}\n");
//	fclose(fp);



    LGM_ARRAY_3D_FREE( Mesh );
    LGM_ARRAY_3D_FREE( Normals );

}
