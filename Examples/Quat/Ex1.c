#include <stdio.h>
#include <Lgm_Quat.h>
#include <Lgm_Vec.h>


int main(void) {

    Lgm_Vector v1, v2, vtmp;
    double q1[4], q2[4], q3[4]; // quaternions

    printf("Start with a vector along x\n");
    Lgm_SetVecElements(&v1, 2, 0, 0);
    Lgm_PrintVector(&v1);
    printf("Rotate it 90 degrees about y\n");
    Lgm_SetVecElements(&vtmp, 0, 1, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q1);
    Lgm_PrintQuat((double *)&q1);
    printf("To get:\n");
    Lgm_QuatRotateVector((double *)&q1, &v1, &v2);
    Lgm_PrintVector(&v2);
    printf("\n\n");

    printf("Start with a vector along x\n");
    Lgm_SetVecElements(&v1, 2, 0, 0);
    Lgm_PrintVector(&v1);
    printf("Rotate it 90 degrees about y\n");
    Lgm_SetVecElements(&vtmp, 0, 1, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q1);
    Lgm_PrintQuat((double *)&q1);
    printf("Rotate it 90 degrees about x\n");
    Lgm_SetVecElements(&vtmp, 1, 0, 0);
    Lgm_AxisAngleToQuat(&vtmp, 90, q2);
    Lgm_PrintQuat((double *)&q2);
    printf("Combine the quaternions\n");
    Lgm_QuatCombineQuats((double *)&q1, (double *)&q2, (double *)&q3);
    printf("To get the quat:\n");
    Lgm_PrintQuat((double *)&q3);
    printf("Which gives the vector:\n");
    Lgm_QuatRotateVector((double *)&q3, &v1, &v2); 
    Lgm_PrintVector(&v2); 

    return 0;

}

