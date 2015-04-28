/* simplex.c */

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "assert.h"
#include "math.h"
#include "simplex2.h"


void swap(double *An, int m, int n, int c1, double *XB, int c2)
{
    //double *temp = (double*)malloc(sizeof(double)*m);
    double temp;

    int idx;

    for (idx = 1; idx <= m; idx++)
    {
        //temp[idx-1] = An[getIndex(m,n, idx, c1)];
        temp = An[getIndex(m,n, idx, c1)];
        An[getIndex(m,n, idx, c1)] = XB[getIndex(m,m, idx, c2)];
        XB[getIndex(m,m, idx, c2)] = temp;
    }


    //free(temp);
}

void printLine(int nFactor)
{
    int n = 0;

    printf("\n");
    for (n = 0; n < nFactor+1; n++)
    {
        printf("-----------------");
    }
    printf("\n");
}


void printMatrix(double *A, int m, int n)
{
    int rowC = 0;
    int colC = 0;

    printf("\nMatrix:\n");

    printLine(n);
    while (rowC < m)
    {
        while (colC < n)
        {
            printf("\t%f", A[getIndex(m, n, rowC+1, colC+1)]);
            colC++;
        }
        //printf(" | %f", b[rowC]);
        printf("\n");
        colC = 0;
        rowC++;
    }

    colC = 0;

    while (colC < n)
    {
        printf("\t");
//        printf("%f", c[colC]);
        colC++;
    }
    //printf(" | %f", zb);

    printf("\n");
    //printf("---------------------------------------\n");
    printLine(n);
}

void printMatrix2(double *An, int m, int n, double *XB, double *Xb, double *c, double *cb, double t, int e, int l)
{
    int rowC = 0;
    int colC = 0;

    printf("\nMatrix (An):\n");

    printLine(n);
    while (rowC < m)
    {
        while (colC < n)
        {
            printf("\t%f", An[getIndex(m, n, rowC+1, colC+1)]);
            colC++;
        }
        printf(" | %f", Xb[rowC]);
        printf("\n");
        colC = 0;
        rowC++;
    }

    rowC = 0;
    colC = 0;

    printf("\nMatrix (XB):\n");

    printLine(n);
    while (rowC < m)
    {
        while (colC < m)
        {
            printf("\t%f", XB[getIndex(m, m, rowC+1, colC+1)]);
            colC++;
        }
//        printf(" | %f", XB[rowC]);
        printf("\n");
        colC = 0;
        rowC++;
    }

    colC = 0;
    printf("c = \n");
    while (colC < n)
    {
        printf("\t");
        printf("%f", c[colC]);
        colC++;
    }

    colC = 0;
    printf("\ncb = \n");
    while (colC < n)
    {
        printf("\t");
        printf("%f", cb[colC]);
        colC++;
    }
    printf("\n");

    printf("t = %f\n",t);
    printf("e = %d\n",e);
    printf("l = %d\n",l);


    //printf("---------------------------------------\n");
    printLine(n);
}

void subtractFromColumn(double *An, int m, int n, int e, const double* t, double *XB, double *Xb)
{
    assert(An);
    assert(XB);
    assert(Xb);
    assert(t);

    int idxBasis = 0;
    int idx = 1;

    for (idxBasis = 0; idxBasis < m; idxBasis++)
    {
        for (idx = 1; idx <= m; idx++)
        {
            double addend = (*t) * An[getIndex(m,n,idx,e)];
            //printf ("Subtracting from basis t(An)e = %f\n A[%d,%d]=%f\t", addend, idx, idxBasis+1, An[getIndex(m,n,idx,e)]);
            //printf("XB{%d,%d} = %f\n", idx, idxBasis+1, XB[getIndex(m, m, idx, idxBasis+1)]);
            XB[getIndex(m, m, idx, idxBasis+1)] = XB[getIndex(m, m, idx, idxBasis+1)] - addend;
            //printf("now XB{%d,%d} = %f\n", idx, idxBasis+1, XB[getIndex(m, m, idx, idxBasis+1)]);
        }
    }
}

void addToColumn(double *An, int m, int n, int e, const double *t)
{
    assert(An);
    assert(e >=1);

    int idx = 1;

    for (idx = 1; idx <= m; idx++)
    {
        //printf ("Updating {%d,%d} current value = %f\n", idx, e, An[getIndex(m,n,idx, e)]);
        An[getIndex(m,n,idx, e)] = An[getIndex(m,n,idx, e)] + (*t);
        //printf ("Updated {%d,%d} value = %f\n", idx, e, An[getIndex(m,n,idx, e)]);
    }
}

void copyRow(double *r, const double* An, int m, int n, int l)
{
    assert(r);
    assert(An);

    int idx;

    for (idx = 1; idx <= n; idx++)
    {
        r[idx-1] = An[getIndex(m,n, l, idx)];
    }
}


void copyColumn(double *d, const double* An, int m, int n, int e, double a)
{
    assert(d);
    assert(An);

    int idx;

    for (idx = 1; idx <= m; idx++)
    {
        d[idx-1] = An[getIndex(m,n, idx, e)]+a;
    }
}

// m and n are are indexes that start at 1
// j is the column to use
// ∥(1,2,3)∥ = 1^2+2^2+3^2 = √14
double calcEuclideanDistance(const double *An, int m, int n, int j)  // line # 2
{
    assert(m > 0 && n > 0);

    int idx = 0;
    double norm = 0;

    for (idx = 0; idx < m; idx++)
    {
        norm += powf(An[getIndex(m, n, idx+1, j)],2);
    }

    norm = sqrtf(norm);

    return norm;
}

// returns the column index (starting at column=1)
// c are the objective function coefficients
int findEnteringVariable(const double *An, int m, int n, const double *c)
{
    assert(An);
    assert(m > 0 && n > 0);

    double max = -1.0;
    double ratio = 0;
    int e = -1;
    int j = 1;

    for (j = 1; j <= n; j++)  // loop through columns and find largest ratio of the c/norm
    {
        double Cj = c[j-1];

        double Gj = calcEuclideanDistance(An, m, n, j);  // line #2

        printf("euclidean distance for col:%d is: %f\n", j, Gj);

        Gj *= Gj; // line #2

        // argmax
        ratio = c[j-1] / sqrtf(Gj); // line #3
        printf("Ratio of {c[%d] = %f} (%f)/(%f) = %f\n",j, c[j-1], c[j-1], sqrtf(Gj), ratio);
        if (ratio >= max)   // line #3
        {
            max = c[j-1]/sqrtf(Gj);
            e = j;  // e stores the column index
        }
    }

    if (max < 0)
    {
        printf ("\n***Optimal Found!\n\n");  // line #5
        exit(0);
    }

    return e;  // // line #3
}


// An (not base matrix
// B  Ax = B
// m rows in
int findLeavingVariable(const double *b, const double *An, int m, int n, int e, int *l, double *t)
{
    *l = -1;
    int rowIdx = 1;
    int idx = 0;
    double new_t = 0;

    for (rowIdx = 1; rowIdx <= m; rowIdx++)
    {
        idx = getIndex(m, n, rowIdx, e);

        if (An[idx] == 0)
        {
            printf ("\n{%d,%d}, has a 0 divisor skipping as unbounded\n", rowIdx, e);
            continue;
        }

        if (*l == -1)
        {
            *t = b[rowIdx-1]/An[idx]; // t= min { bi/Aie | cj > 0 for all j element of N
            *l = rowIdx;
            printf ("{%d,%d}, is new min with t = %f\n", rowIdx, e, *t);
        }
        else
        {
            new_t = b[rowIdx-1]/An[idx];

            if (new_t < *t)
            {
                *t = new_t;
                *l = rowIdx;
            }
        }
    }

    if ( *l == -1 )
    {
        printf ("\n***Unbounded!\n");
        exit(0);
    }
}

int main (int argc, char *argv[])
{
    int m = 3;  // # of rows and          size of columns
    int n = 2;  // # of cols and          size of rows
    int e = -1; // entering column
    int l = -1; // leaving row
    int lda = 3; /* Leading dimension of 5 * 4 matrix is 5 */

    double *An = (double*)malloc(sizeof(double)*m*n);
    double *XB = (double*)malloc(sizeof(double)*m*m);

    /*
     An =
      2 3
      1 1
      0 1

     c = 2 3
l
     b = 20
          4
         10

     XB =
        1 0 0
        0 1 0
        0 0 1

     Xb =
          20
           4
           3

     cb = 0 0
     */

    /* The elements of the first column */
    An[0] = 2;
    An[1] = 1;
    An[2] = 0;

    XB[0] = 1;
    XB[1] = 0;
    XB[2] = 0;
    //a[3] = 4;
    /* The elements of the second column */
    An[m] = 3;
    An[m+1] = 1;
    An[m+2] = 1;
    //a[m+3] = 1;

    XB[m] = 0;
    XB[m+1] = 1;
    XB[m+2] = 0;

    /* The elements of the third column */
    //An[m*2] = 0;
    //An[m*2+1] = 1;
    //a[m*2+2] = 0;
    //a[m*2+3] = 6;

    XB[m*2] = 0;
    XB[m*2+1] = 0;
    XB[m*2+2] = 1;

    /* The elements of the fourth column */
    //a[m*3] = 0;
    //a[m*3+1] = 1;
    //a[m*3+2] = 0;
    //a[m*3+3] = 8;

    double *c = (double*)malloc(sizeof(double)*n);
    c[0] = 2;
    c[1] = 4;
    /*
    c[3] = 0;*/

    double *cB = (double*)malloc(sizeof(double)*m);
    cB[0] = 0;
    cB[1] = 0;
    cB[2] = 0;

    double *b = (double*)malloc(sizeof(double)*m);
    b[0] = 20;
    b[1] = 4;
    b[2] = 10;

    // set these variables
    double *Xb = (double*)(malloc(sizeof(double)*m));   // slack variables and negative offsets are in the basis initially. values in array are column indexes
    Xb[0] = 20;
    Xb[1] = 4;
    Xb[2] = 3;


    printMatrix(An,m,n);
    printMatrix(XB,m,m);

    int iter = 0;
    while(iter < 100)
    {
        e = findEnteringVariable(An, m, n, c);

        printf("Entering column = %d\n", e);

        double t = -1.0;
        findLeavingVariable(b, An, m, n, e, &l, &t);
        //printf("Leaving row = %d\n", l);
        //printf("Pivot value = %f {%d,%d}\n", An[getIndex(m,n,l,e)], l, e);
        //printf("t = %f\n", t);

        double *d = (double*)malloc(sizeof(double)*m);

        copyColumn(d, An, m, n, e, 0); // line # 13
        d[e-1] = An[getIndex(m,n,l,e)]-1;

        printf("d**T: ->\n");
        printArray(1,m,d,1);

        double *r = (double*)malloc(sizeof(double)*n);

        printf("r: ->\n");
        copyRow(r, An, m, n, l); // line # 14

        printArray(1,n,r,1);

        printf("\nt = %f\n",t);
        addToColumn(An, m, n, e,  &t); // line #15

        printf("\nLine #15\n");
        printMatrix(An,m,n);

        printf("\nLine #16\n");
        subtractFromColumn(An, m, n, e, &t, XB, Xb); // line #16

        printMatrix(XB,m,m);


        enum CBLAS_ORDER order;
        order = CblasColMajor;

        double alpha = -1.0/An[getIndex(m,n,l,e)];
        printf ("alpha = %f\n",alpha);
        int incx = 1;
        int incy = 1;

        printf("\nbefore dfger:->");
        printMatrix(An,m,n);
        printf("d**T: ->\n");
        printArray(1,m,d,1);
        printf("r: ->\n");
        printArray(1,n,r,0);


        cblas_dger(order, m, n, alpha, r, incx, d, incy, An, lda);
        printf("\nafter dfger:->");
        //printMatrix(An, m, n);
        printMatrix2(An,m,n,XB,Xb,c,cB,t,e,l);

        alpha = -c[l-1]/An[getIndex(m,n,l,e)];
        printf ("alpha = %f\n",alpha);
        cblas_daxpy(n, alpha, r, incx, c, incy);
        printf("\nafter daxby:->");
        //printMatrix(An, m, n);
        printMatrix2(An,m,n,XB,Xb,c,cB,t,e,l);

        swap(An, m, n, e, XB, l);

        double temp = c[e-1];
        c[e-1] = cB[l-1];
        cB[l-1] = temp;

        printf("\nafter swap:->");
        //printMatrix(An, m, n);
        printMatrix2(An,m,n,XB,Xb,c,cB,t,e,l);




        iter++;
    }
    return 0;
}

