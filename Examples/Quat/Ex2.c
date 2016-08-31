#include <stdio.h>
#include <Lgm_Quat.h>
#include <Lgm_Vec.h>


int main(void) {

    double matrix[3][3];
    double q1[4], q2[4], q3[4];
    Lgm_Vector v1, v2, v3, v4, vtmp;
    
    printf("Start with a vector along x\n");
    Lgm_SetVecElements(&v1, 2, 0, 0);
    Lgm_SetVecElements(&v3, 2, 0, 0); // used later
    Lgm_PrintVector(&v1);

    printf("Make Quat that rotates 90 deg about y\n");
    Lgm_SetVecElements(&vtmp, 0, 1, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q1);
    Lgm_PrintQuat((double *)&q1);
    printf("To get:\n");
    Lgm_QuatRotateVector((double *)&q1, &v1, &v2);
    Lgm_PrintVector(&v2);
    printf("\nMake that quat into a rotation matrix\n");
    Lgm_Quat_To_Matrix(q1, matrix);
    printf("%lf %lf %lf\n", matrix[0][0],matrix[0][1],matrix[0][2]);
    printf("%lf %lf %lf\n", matrix[1][0],matrix[1][1],matrix[1][2]);
    printf("%lf %lf %lf\n", matrix[2][0],matrix[2][1],matrix[2][2]);
    printf("And use this matrix to rotate the vector, to get:\n");
    Lgm_MatTimesVec(matrix, &v3, &v4);
    Lgm_PrintVector(&v4);

    printf("\n\n\n\n");

    printf("Start with a vector along x\n");
    Lgm_SetVecElements(&v1, 2, 0, 0);
    Lgm_SetVecElements(&v3, 2, 0, 0); // used later
    Lgm_PrintVector(&v1);

    printf("Make Quat that rotates 90 deg about y\n");
    Lgm_SetVecElements(&vtmp, 0, 1, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q1);
    Lgm_PrintQuat((double *)&q1);
    printf("And a Quat that rotates 90 deg about x\n");
    Lgm_SetVecElements(&vtmp, 1, 0, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q2);
    Lgm_PrintQuat((double *)&q2);
    printf("Combine them to get a single rotation\n");
    Lgm_QuatCombineQuats(q1, q2, q3);
    Lgm_PrintQuat((double *)&q3);

    printf("Use this on the vector to get:\n");
    Lgm_QuatRotateVector((double *)&q3, &v1, &v2);
    Lgm_PrintVector(&v2);
    printf("\nMake that quat into a rotation matrix\n");
    Lgm_Quat_To_Matrix(q3, matrix);
    printf("%lf %lf %lf\n", matrix[0][0],matrix[0][1],matrix[0][2]);
    printf("%lf %lf %lf\n", matrix[1][0],matrix[1][1],matrix[1][2]);
    printf("%lf %lf %lf\n", matrix[2][0],matrix[2][1],matrix[2][2]);
    printf("And use this matrix to rotate the vector, to get:\n");
    Lgm_MatTimesVec(matrix, &v3, &v4);
    Lgm_PrintVector(&v4);

    printf("\nOr Start with the rotation matrix and make it a Quat\n");
    Lgm_MatrixToQuat(matrix, q3);
    printf("Which gives:\n");
    Lgm_PrintQuat((double *)&q3);
   


    

    
    return 0;

}

