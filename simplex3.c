#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "assert.h"
#include "string.h"
#include "math.h"

int getIndex(int nRows, int nCols, int r, int c)
{
    r--;
    c--;

    assert(r>=0 && r < nRows);
    assert(c>=0 && c < nCols);
    return c*nRows+r;
}

void swapColumn(double *An, int m, int n, int c1, double *XB, int c2)
{
    double temp;

    int idx;

    for (idx = 1; idx <= m; idx++)
    {
        temp = An[getIndex(m,n, idx, c1)];
        An[getIndex(m,n, idx, c1)] = XB[getIndex(m,n, idx, c2)];
        XB[getIndex(m,n, idx, c2)] = temp;
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

void copyRow2(const double *r, double* An, int m, int n, int l)
{
    assert(r);
    assert(An);

    int idx;

    for (idx = 1; idx <= n; idx++)
    {
        An[getIndex(m,n, l, idx)] = r[idx-1];
    }
}


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

void eliminate(double *a, int m, int n, int e, int l, double *b, double *c, double *z)
{
    enum CBLAS_ORDER order;

    int lda = m; /* Leading dimension of 5 * 4 matrix is 5 */
    int incx = 1;
    int incy = 1;

    order = CblasColMajor;
    double alpha = -1.0 * ((a[getIndex(m,n,l,e)]-1)/(a[getIndex(m,n,l,e)]));

    double *x = (double *)malloc(sizeof(double)*m);
    double *y = (double *)malloc(sizeof(double)*n);
    double *orig_y = (double *)malloc(sizeof(double)*n);
    double *orig_b = (double *)malloc(sizeof(double)*m);
    double *temp_r = (double *)malloc(sizeof(double)*n);
    double *temp_c = (double *)malloc(sizeof(double)*m);

    memset(x,0, sizeof(double)*m);
    int idx,idx2;

    for (idx = 1; idx <= n; idx++)
    {
        y[idx-1] = a[getIndex(m,n,l,idx)];
        orig_y[idx-1] = a[getIndex(m,n,l,idx)];
    }

    for (idx = 1; idx <= m; idx++)
    {
       orig_b[idx-1] = b[idx-1];
    }

    x[l-1] = 1;

    /*printf("1st:\n");
    printArray(1,n,y,0,"y");
    printArray(1,m,x,0,"x");
    printArray(1,m,b,0,"b");
    printArray(1,m,orig_b,0,"orig_b");*/
    cblas_dger(order, m, n, alpha, x, incx, y, incy, a, lda);
    b[l-1] += alpha*b[l-1];

    printMatrix(a, m, n, b, c, *z, "set pivot to 1");
    /***
     * zero out the rest of equations
     * */
    for (idx = 1; idx <= m; idx++)
    {
        if (idx == l)
            continue;

        for (idx2 = 1; idx2 <= n; idx2++)
        {
            y[idx-1] = a[getIndex(m,n,idx,idx2)];
        }

        memset(x,0, sizeof(double)*m);
        x[idx-1] = 1;

        alpha = -1.0 / ( orig_y[e-1]/(a[getIndex(m,n,idx, e)]));
        printf("alpha = %f\n",alpha);

        b[idx-1] += (alpha * orig_b[l-1]);

        /*printf("2nd:\n");
        printf("alpha = %f\n", alpha);
        printArray(1,n,orig_y,0,"orig_y");
        printArray(1,n,y,0,"y");
        printArray(1,m,x,0,"x");
        printArray(1,m,b,0,"b");
        printArray(1,m,orig_b,0,"orig_b");*/

        cblas_dger(order, m, n, alpha, x, incx, orig_y, incy, a, lda);
        //printArray(1,m,b,0,"b");
    }

    double multiple = -1.0 * ((c[e-1]/(a[getIndex(m,n,l,e)])));

    copyRow(temp_r, a, m, n, l);

    printf("multiple = %f\n", multiple);
    /*printf("temp_r = \n");
    printArray(1, n, temp_r, 0,"c");
    printf("c = \n");
    printArray(1, n, c, 0,"c");*/

    cblas_daxpy(n, multiple, temp_r, incx, c, incy);

    /*printf("temp_r = \n");
    printArray(1, n, temp_r,0,"temp_r");
    printf("c = \n");
    printArray(1, n, c, 0,"c");*/

    //printf("b[l-1] = %f\n", b[l-1]);

/*    printf("b =\n");
    printArray(m,1, b, 1, "b");*/
    *z = *z + multiple * b[l-1];
    printf ("z = %f\n", *z);

    free(x);
    free(y);
    free(temp_r);
    free(temp_c);
}

void printArray(int m, int n, double *c, int transpose, const char* label)
{
    int rowC = 0;
    int colC = 0;

    while (rowC < m)
    {
        while(colC < n)
        {
            if (transpose == 0)
            {
                printf (" %s[%d,%d] = %f\t", label,rowC+1, colC+1, c[getIndex(m, n, rowC+1, colC+1)]);
            }
            else
            {
                printf (" %s[%d,%d] = %f\n", label,rowC+1, colC+1, c[getIndex(m, n, rowC+1, colC+1)]);
            }
            colC++;
        }
        printf("\n");
        colC=0;
        rowC++;
    }
}

void printLine(int nFactor)
{
    int n = 0;

    for (n = 0; n < nFactor+1; n++)
    {
        printf("-----------------");
    }
    printf("\n");
}


void printMatrix(double *a, int m, int n, double *b, double *c, double z, const char *label)
{
    int rowC = 0;
    int colC = 0;

    printf("\nMatrix (%s):\n",label);

    printLine(n);
    while (rowC < m)
    {
        while (colC < n)
        {
            if (a[getIndex(m, n, rowC+1, colC+1)] >= 0)
            {
                printf("\t %f",a[getIndex(m, n, rowC+1, colC+1)]);
            }
            else
            {
                printf("\t%f", a[getIndex(m, n, rowC+1, colC+1)]);
            }
            colC++;
        }
        printf ("\t|    %f\n",b[rowC]);

        if (rowC+1 < m)
        {
         //   printf("\n");
        }
        colC = 0;
        rowC++;
    }
    printLine(n);

    for (colC = 1; colC <= n; colC++)
    {
        printf("\t%f", c[colC-1]);
    }
    printf("\t|   %f",z);
    printf("\n");
    printLine(n);
    printf("\n");
}

// returns the column index (starting at column=1)
// c are the objective function coefficients
int findEnteringVariable(const double *a, int m, int n, const double *c)
{
    assert(a);
    assert(m > 0 && n > 0);

    double max = -1.0;
    double ratio = 0;
    int first = 1;
    int e = -1;
    int j = 1;

    for (j = 1; j <= n; j++)  // loop through columns and find largest ratio of the c/norm
    {
        double Cj = c[j-1];

        double Gj = calcEuclideanDistance(a, m, n, j);  // line #2

        //printf("euclidean distance for col:%d is: %f\n", j, Gj);

        Gj *= Gj; // line #2

        // argmax
        ratio = c[j-1] / sqrtf(Gj); // line #3
        //printf("Ratio of {c[%d] = %f} (%f)/(%f) = %f\n",j, c[j-1], c[j-1], sqrtf(Gj), ratio);
        if (ratio >= max || first == 1)   // line #3
        {
            max = c[j-1]/sqrtf(Gj);
            e = j;  // e stores the column index
            first = 0;
        }
    }

    if (max > 0)
    {
        printf ("\n***Optimal Found!\n\n");  // line #5
        exit(0);
    }

    printf ("Entering column = %d\n",e);
    return e;  // // line #3
}


int findEnteringVariable2(const double *a, int m, int n, const double *c)
{
    assert(a);
    assert(m > 0 && n > 0);

    double max = -1.0;
    double ratio = 0;
    int allpos = 1;
    int first = 1;
    int e = -1;
    int j = 1;

    for (j = 1; j <= n; j++)  // loop through columns and find largest ratio of the c/norm
    {
        if (c[j-1] < 0)
        {
            allpos = -1;
        }
        ratio = abs(c[j-1]);
        if (ratio >= max || first == 1)   // line #3
        {
            max = abs(ratio);
            e = j;  // e stores the column index
            first = 0;
        }
    }

    if (allpos > 0)
    {
        printf ("\n***Optimal Found!\n\n");  // line #5
        exit(0);
    }

    printf ("Entering column = %d\n",e);
    return e;  // // line #3
}

int findEnteringVariable3(const double *a, int m, int n, const double *c)
{
    assert(a);
    assert(m > 0 && n > 0);

    double min = -1.0;
    double ratio = 0;
    int allpos = 1;
    int first = 1;
    int e = -1;
    int j = 1;

    for (j = 1; j <= n; j++)  // loop through columns and find largest ratio of the c/norm
    {
        if (c[j-1] < 0)
        {
            allpos = -1;
        }
        ratio = c[j-1];
        if (ratio <= min || first == 1)   // line #3
        {
            min = ratio;
            e = j;  // e stores the column index
            first = 0;
        }
    }

    if (allpos > 0)
    {
        printf ("\n***Optimal Found!\n\n");  // line #5
        exit(0);
    }

    printf ("Entering column = %d\n",e);
    return e;  // // line #3
}



// An (not base matrix
// B  Ax = B
// m rows in
int findLeavingVariable(const double *b, const double *a, int m, int n, int e, int *l, double *t)
{
    *l = -1;
    int rowIdx = 1;
    int idx = 0;
    double new_t = 0;

    for (rowIdx = 1; rowIdx <= m; rowIdx++)
    {
        idx = getIndex(m, n, rowIdx, e);

        if (a[idx] == 0)
        {
            printf ("\n{%d,%d}, has a 0 divisor skipping as unbounded\n", rowIdx, e);
            continue;
        }
        if (a[idx] < 0)
        {
            printf ("\n{%d,%d}, is < 0 skipping\n", rowIdx, e);
            continue;
        }

        if (*l == -1)
        {
            *t = b[rowIdx-1]/a[idx]; // t= min { bi/Aie | cj > 0 for all j element of N
            *l = rowIdx;
            printf ("{%d,%d}, is new min with t = %f\n", rowIdx, e, *t);
        }
        else
        {
            new_t = b[rowIdx-1]/a[idx];

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

    printf("Leaving row = %d\n", *l);
}


int main ( )
{
   enum CBLAS_ORDER order;
   enum CBLAS_TRANSPOSE transa;

   double *a, *b, *c;
   double t;
   double alpha, beta;
   double z = 0;
   int m, n, lda, incx, incy, i;

   order = CblasColMajor;
   transa = CblasNoTrans;

   m = 3; /* Size of Column ( the number of rows ) */
   n = 5; /* Size of Row ( the number of columns ) */
   lda = 3; /* Leading dimension of 5 * 4 matrix is 5 */
   incx = 1;
   incy = 1;
   alpha = 1;
   beta = 0;

   a = (double *)malloc(sizeof(double)*m*(n+m));
   b = (double *)malloc(sizeof(double)*n);
   c = (double *)malloc(sizeof(double)*(n+m));

   b[0] = 66;
   b[1] = 180;
   b[2] = 200;
   //b[3] = 1;

   c[0] = -90;
   c[1] = -75;
   c[2] = 0;
   c[3] = 0;
   c[4] = 0;
   //c[5] = 0;

   /* The elements of the first column */
   a[0] = 3;
   a[1] = 9;
   a[2] = 2;
   //a[3] = 4;
   /* The elements of the second column */
   a[m] = 2;
   a[m+1] = 4;
   a[m+2] = 10;
   //a[m+3] = 8;
   /* The elements of the third column */
   a[m*2] = 1;
   a[m*2+1] = 0;
   a[m*2+2] = 0;
   //a[m*2+3] = 12;
   /* The elements of the fourth column */
   a[m*3] = 0;
   a[m*3+1] = 1;
   a[m*3+2] = 0;
   //a[m*3+3] = 16;

   /* The elements of the fifth column */
   a[m*4] = 0;
   a[m*4+1] = 0;
   a[m*4+2] = 1;

   int iter = 0;
   while (iter < 50)
   {
       printMatrix(a, m, n, b, c, z, "initial matrix");
       //printArray(1, n, c, 0,"c");

       //int e = findEnteringVariable2(a, m, n-m, c);
       int e = findEnteringVariable3(a, m, n, c);
       int l = -1;
       //findLeavingVariable(b, a, m ,n-m, e, &l, &t);
       findLeavingVariable(b, a, m ,n, e, &l, &t);

       eliminate(a, m, n, e, l, b, c, &z);

       //int val = m+l;
       //printf("val =%d\n",val);
       //swapColumn(a, m, n, e, a, val);
       //swapColumn(a, 1, n, e, c, val);

       printMatrix(a, m, n, b, c, z, "eliminated matrix");
       //swapColumn(a,m,n,e,a,m+l);
       //swapColumn(c,1,n,e,c,m+l);
       //printMatrix(a, m, n, b, c, z, "swapped matrix");
       //printArray(1, n, c, 0,"c");
       iter++;
    }

//   eliminate(a, m, n, 1, 1);
//   printMatrix(a, m, n);

   free(a);
   free(b);
   free(c);
   return 1;
}
