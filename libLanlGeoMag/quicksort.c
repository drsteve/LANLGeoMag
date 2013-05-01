#include <stdio.h>
#include <stdlib.h>
#include "Lgm/quicksort.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M	7
#define NSTACK 2001

/*
 * These routines assume the retarded NR index starts at 1 thing.
 */
void quicksort( unsigned long n, double *arr ) {


    unsigned long 	i, ir=n, j, k, l=1;
    int 		jstack=0, istack[NSTACK];
    double		a, temp;


    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++){
                a = arr[j];
                for (i=j-1; i>=1; i--){
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                }
                arr[i+1] = a;
            }
            if (jstack == 0) break;
            ir = istack[jstack--];
            l  = istack[jstack--];
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            if ( arr[l+1] > arr[ir] ) {
                SWAP(arr[l+1],arr[ir]);
            }
            if ( arr[l] > arr[ir] ) {
                SWAP(arr[l],arr[ir]);
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l]);
            }
            i = l+1;
            j = ir;
            a = arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
            }
            arr[l] = arr[j];
            arr[j] = a;
            jstack += 2;
            if ( jstack > NSTACK-1 ) {
                printf("NSTACK too small in quicksort\n");
                exit(-1);
            }
            if ( ir-i+1 >= j-1 ) {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack]   = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }



}


void quicksort_uli( unsigned long n, unsigned long *arr ) {


    unsigned long 	i, ir=n, j, k, l=1;
    int 		    jstack=0, istack[NSTACK];
    unsigned long	a, temp;


    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++){
                a = arr[j];
                for (i=j-1; i>=1; i--){
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                }
                arr[i+1] = a;
            }
            if (jstack == 0) break;
            ir = istack[jstack--];
            l  = istack[jstack--];
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            if ( arr[l+1] > arr[ir] ) {
                SWAP(arr[l+1],arr[ir]);
            }
            if ( arr[l] > arr[ir] ) {
                SWAP(arr[l],arr[ir]);
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l]);
            }
            i = l+1;
            j = ir;
            a = arr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
            }
            arr[l] = arr[j];
            arr[j] = a;
            jstack += 2;
            if ( jstack > NSTACK-1 ) {
                printf("NSTACK too small in quicksort\n");
                exit(-1);
            }
            if ( ir-i+1 >= j-1 ) {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack]   = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }



}



void quicksort2(  unsigned long n, double *arr, double *brr ) {


    unsigned long 	i, ir=n, j, k, l=1;
    int 		jstack=0, istack[NSTACK];
    double		a, b, temp;


    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++){
                a = arr[j];
                b = brr[j];
                for (i=j-1; i>=1; i--){
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
            }
            if (!jstack) {
                return;
            }
            ir = istack[jstack];
            l  = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            SWAP(brr[k],brr[l+1]);
            if ( arr[l+1] > arr[ir] ) {
                SWAP(arr[l+1],arr[ir]);
                SWAP(brr[l+1],brr[ir]);
            }
            if ( arr[l] > arr[ir] ) {
                SWAP(arr[l],arr[ir]);
                SWAP(brr[l],brr[ir]);
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l]);
                SWAP(brr[l+1],brr[l]);
            }
            i = l+1;
            j = ir;
            a = arr[l];
            b = brr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
                SWAP(brr[i],brr[j]);
            }
            arr[l] = arr[j];
            arr[j] = a;
            brr[l] = brr[j];
            brr[j] = b;
            jstack += 2;
            if ( jstack > NSTACK-1 ) {
                printf("NSTACK too small in quicksort2\n");
                exit(-1);
            }
            if ( ir-i+1 >= j-1 ) {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack]   = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}

void quicksort2uli(  unsigned long n, double *arr, unsigned long int *brr ) {


    unsigned long 	i, ir=n, j, k, l=1;
    int 		jstack=0, istack[NSTACK];
    int         b;
    double		a, temp;


    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++){
                a = arr[j];
                b = brr[j];
                for (i=j-1; i>=1; i--){
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
            }
            if (!jstack) {
                return;
            }
            ir = istack[jstack];
            l  = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            SWAP(brr[k],brr[l+1]);
            if ( arr[l+1] > arr[ir] ) {
                SWAP(arr[l+1],arr[ir]);
                SWAP(brr[l+1],brr[ir]);
            }
            if ( arr[l] > arr[ir] ) {
                SWAP(arr[l],arr[ir]);
                SWAP(brr[l],brr[ir]);
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l]);
                SWAP(brr[l+1],brr[l]);
            }
            i = l+1;
            j = ir;
            a = arr[l];
            b = brr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
                SWAP(brr[i],brr[j]);
            }
            arr[l] = arr[j];
            arr[j] = a;
            brr[l] = brr[j];
            brr[j] = b;
            jstack += 2;
            if ( jstack > NSTACK-1 ) {
                printf("NSTACK too small in quicksort2i\n");
                exit(-1);
            }
            if ( ir-i+1 >= j-1 ) {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack]   = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}

void quicksort3(  unsigned long n, double *arr, double *brr, double *crr ) {


    unsigned long 	i, ir=n, j, k, l=1;
    int 		    jstack=0, istack[NSTACK];
    double		    a, b, c, temp;


    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++){
                a = arr[j];
                b = brr[j];
                c = crr[j];
                for (i=j-1; i>=1; i--){
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                    crr[i+1] = crr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
                crr[i+1] = c;
            }
            if (!jstack) {
                return;
            }
            ir = istack[jstack];
            l  = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k],arr[l+1]);
            SWAP(brr[k],brr[l+1]);
            SWAP(crr[k],crr[l+1]);
            if ( arr[l+1] > arr[ir] ) {
                SWAP(arr[l+1],arr[ir]);
                SWAP(brr[l+1],brr[ir]);
                SWAP(crr[l+1],crr[ir]);
            }
            if ( arr[l] > arr[ir] ) {
                SWAP(arr[l],arr[ir]);
                SWAP(brr[l],brr[ir]);
                SWAP(crr[l],crr[ir]);
            }
            if (arr[l+1] > arr[l]) {
                SWAP(arr[l+1],arr[l]);
                SWAP(brr[l+1],brr[l]);
                SWAP(crr[l+1],crr[l]);
            }
            i = l+1;
            j = ir;
            a = arr[l];
            b = brr[l];
            c = crr[l];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
                SWAP(brr[i],brr[j]);
                SWAP(crr[i],crr[j]);
            }
            arr[l] = arr[j];
            arr[j] = a;
            brr[l] = brr[j];
            brr[j] = b;
            crr[l] = crr[j];
            crr[j] = c;
            jstack += 2;
            if ( jstack > NSTACK-1 ) {
                printf("NSTACK too small in quicksort2\n");
                exit(-1);
            }
            if ( ir-i+1 >= j-1 ) {
                istack[jstack]   = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack]   = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}


/*
 * An implementation of bubbleSort in case a user wants bubbleSort instead of quicksort
 * BAL - 20Apr2011
 */
void bubbleSort(unsigned long n, double *arr)
{
  long i, j, temp;

  for (i=(n-1); i>0; i--) /* loop over first element */
  {
    for (j=1; j<=i; j++) /* loop over second element */
    {
      if (arr[j-1] > arr[j])
      {
        SWAP(arr[j-1],arr[j])
      }
    }
  }
}
