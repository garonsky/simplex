/* simplex.c */

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "assert.h"
#include "math.h"

int foundOptimal = 0;


//int getEnteringVariable(const double An,


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

// m : rows
// n : cols
// a :
void printMatrix(int m, int n, double *a, double *b, double *c, double zb)
{
    int rowC = 0;
    int colC = 0;

    printf("\nMatrix:\n");

    printLine(n);
    while (rowC < m)
    {
        while (colC < n)
        {
            printf("\t%f", a[getIndex(m, n, rowC+1, colC+1)]);
            colC++;
        }
        printf(" | %f", b[rowC]);
        printf("\n");
        colC = 0;
        rowC++;
    }

    colC = 0;

    while (colC < n)
    {
        printf("\t");
        printf("%f", c[colC]);
        colC++;
    }
    printf(" | %f", zb);

    printf("\n");
    //printf("---------------------------------------\n");
    printLine(n);
}


void printArray(int m, int n, double a[])
{
    int rowC = 0;
    int colC = 0;

    while (rowC < m)
    {
        while(colC < n)
        {
            printf (" a[%d,%d] = %f\t", rowC+1, colC+1, a[getIndex(m, n, rowC+1, colC+1)]);
            colC++;
        }
        printf("\n");
        colC=0;
        rowC++;
    }
}

int getIndex(int nRows, int nCols, int r, int c)
{
    r--;
    c--;

    assert(r>=0);
    assert(c>=0);
    return c*nRows+r;
}

int getMaxCol2(double a[], int nRows, int nCols, int nStartRow)
{
    int foundPos= 0;
    int currentIdx = getIndex(nRows,nCols, nStartRow, 1);
    int endIdx = getIndex(nRows, nCols, nStartRow, nCols);

    double max = -1;

    int i, max_i, arrIdx = -1;

    for (i = 1; i <= nCols; i += 1)
    {
        arrIdx = getIndex(nRows, nCols, nStartRow, i);
        if  (a[arrIdx] > 0)
        {
            foundPos = 1;
        }

        if (max < a[arrIdx])
        {
            max = a[arrIdx];
            max_i = i;
        }
    }

    if (foundPos == 0)
    {
        //printf("Optimal found!");
        foundOptimal = 1;
    }

    return max_i;
}

int getMinRow2(double a[], int nRows, int nCols, int startCol, double *b)
{
    int nStartRow = 1;
    int endIdx = getIndex(nRows, nCols, nStartRow, nCols);

    double max = 0;

    int i;
    int max_i = nStartRow;
    int arrIdx = -1;

    for (i = nStartRow; i <= nRows; i += 1)
    {
        arrIdx = getIndex(nRows, nCols, i, startCol);
        endIdx = getIndex(nRows, nCols, nStartRow, nCols);
        if (a[arrIdx] > 0 && max < a[endIdx]/b[i])
        {
            max = a[endIdx]/b[i]);
            max_i = i;
        }
    }

    return max_i;
}

int main (int argc, char *argv[])
{
   enum CBLAS_ORDER order;
   enum CBLAS_TRANSPOSE transa;

   double *a, *x, *y, *b, *c;
   double alpha, beta, zb;
   int m, n, lda, incx, incy, i;

   order = CblasColMajor;
   transa = CblasNoTrans;
   int iteration = 1;

   m = 3; /* Size of Column ( the number of rows ) */
   n = 5; /* Size of Row ( the number of columns ) */
   lda = 5; /* Leading dimension of 5 * 4 matrix is 5 */
   incx = 1;
   incy = 1;
   alpha = 1;
   beta = 0;
   int IdentityMatrixOffSet = 2;


   a = (double *)malloc(sizeof(double)*m*n);
   x = (double *)malloc(sizeof(double)*n);
   y = (double *)malloc(sizeof(double)*n);
   c = (double *)malloc(sizeof(double)*n);
   b = (double *)malloc(sizeof(double)*m);

   zb = 0.0;
   b[0] = 20;
   b[1] = 4;
   b[2] = 10;

   c[0] = 2;
   c[1] = 4;
   c[2] = 0;
   c[3] = 0;
   c[4] = 0;
   c[5] = 0;
   //c[6] = 0;


   /* The elements of the first column */
   a[0] = 2;
   a[1] = 1;
   a[2] = 0;
   //a[3] = 4;
   /* The elements of the second column */
   a[m] = 3;
   a[m+1] = 1;
   a[m+2] = 2;
   //a[m+3] = 1;
   /* The elements of the third column */
   a[m*2] = 1;
   a[m*2+1] = 0;
   a[m*2+2] = 0;
   //a[m*2+3] = 6;
   /* The elements of the fourth column */
   a[m*3] = 0;
   a[m*3+1] = 1;
   a[m*3+2] = 0;
   //a[m*3+3] = 8;

    /* The elements of the fifth column */
   a[m*4] = 0;
   a[m*4+1] = 0;
   a[m*4+2] = 1;
   //a[m*4+3] = 8;

   /*a[m*5] = 0;
   a[m*5+1] = 0;
   a[m*5+2] = 1;*/
   //a[m*5+3] = 8;

   //a[m*6] = 0;
   //a[m*6+1] = 10;
   //a[m*6+2] = 15;
   //a[m*6+3] = 8;


   printf("\nBefore:\n");
   printMatrix(m, n, a, b, c, zb);


   int e = 1;
   int l = 1;

   int iter  = 1;
   if (argc == 1)
   {
       iter = 1;
   }
   else
   {
       iter  = atoi(argv[1]);
   }
   while (iter-- > 0 && foundOptimal == 0)
   {
       printf("\nIteration = %d\n\n",iteration++);
       // Find entering variable
       e = getMaxCol2(c, 1, n, 1);   // e col index

       printf("entering column (e) = %d\t", e);

       // Find leaving variable
       l = getMinRow2(a, m, n, e, b);  // l row index

       printf("leaving row (l) = %d\t", l);

       printf ("pivot value = %f\n", a[getIndex(m,n,l,e)]);

       if (l < 0)
       {
           printf ("\nUnbounded!\n");
           exit(0);
       }

       // t = min ratio of all rows in col e (which is l) from getMinRow
       double t = b[l] / a[getIndex(m,n, l, e)];
       printf("t = %f\n", t);

       int temp = 0;

       // pivoting
       double *d =  (double *)malloc(sizeof(double)*m);
       for (temp = 0; temp < m; temp++)
       {
           d[temp] = a[getIndex(m, n, 1 + temp, e)];
       }
       d[e-1] = a[getIndex(m,n,l,e)] - 1;

       double *r = (double *)malloc(sizeof(double)*(n-IdentityMatrixOffSet));
       for (temp = IdentityMatrixOffSet; temp < n; temp++)
       {
           r[temp-IdentityMatrixOffSet] = a[getIndex(m, n, l, temp+1)];
       }

       // Xe <-- Xe + t
       c[e] += t;
       printf("\nXe (c[e]) = %f\t", c[e]);

       //zb -= (t * c[e]);
       //printf("Xb (zb) = %f\n", zb);
       for (temp = IdentityMatrixOffSet; temp < n; temp++)
       {
           a[getIndex(m, n, m, temp+1)] -= (t * c[e]);
       }


       alpha = 1.0/a[getIndex(m,n,l,e)];

       printf("\nd:\n");
       printArray(1, m, d);
       printf("\nr:\n");
       printArray(1, n-IdentityMatrixOffSet, r);

       printf("\nAlpha = %f\n", alpha);

       cblas_dger(order, m, n, alpha, d, incx, r, incy, a, m);

       alpha *= -c[l];
       printf("Alpha = %f\n", alpha);

       //cblas_daxpy(2, alpha, r, incx, &(a[getIndex(m,n,1,n)]), incy);
       cblas_daxpy(2, alpha, r, incx, c, incy);

       //printf("\nIteration %d:\n", iter);
       double tempE = c[l];
       c[l] = c[e];
       c[e] = tempE;

       printMatrix(m,n,a,b,c,zb);
   }

   if (foundOptimal == 1)
   {
        printf("\n****Optimal found!\n");
   }

   //free(a);
   //free(x);
   //free(y);
   return 1;
}
