/* simplex.c */

#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"
#include "assert.h"
#include "math.h"
#include "simplex2.h"


// m and n are are indexes that start at 1
// c is the column to use
// ∥(1,2,3)∥ = 1^2+2^2+3^2 = √14
double calcEuclideanDistance(const double *d, int m, int n, int c)
{
    assert(d);
    assert(len > 0);

    int idx = 0;
    int norm = 0;

    for (idx = 0; idx < len; idx++)
    {
        norm =+ d[idx] * d[idx];
    }

    norm = sqrtf(norm);

    return norm;

}

// returns the column index (starting at column=1)l
int getEnteringVariable(const double *An, int nRows, int nCols)
{
    assert(An);
    assert(nRows > 0 && nCols > 0);

    double max = -1.0;
    double c;
    int e = -1;
    int j = 1;


    for (j = 1; j <= nCols; j++)
    {
        c = An[getIndex(nRows, nCols, nRows, j)];  // coefficients cAx

        double Gj = calcEuclideanDistance(An, nRows, nCols, j);  // line #2
        Gj *= Gj;

        if (sqrt(Gj) >= max)
        {
            e = j;  // e stores the column index
        }
    }

    if (max < 0)
    {
        printf ("\n***Optimal Found!");
        exit(0);
    }

    return max;  // // line #3
}


// An (not base matrix
// B  Ax = B
// m rows in
void expand(const double *B, const double An, int m, int n, int e, int &l, double &t)
{
    int rowIdx = 1;
    double new_t = 0;

    for (rowIdx = 1; rowIdx <= m; rowIdx++)
    {
        idx = getIndex(m, n, rowIdx, e);

        if (rowIdx == 1)
        {
            t = An[idx]/B[rowIdx--];
            l = rowIdx; // 1
        }
        else
        {
            new_t = An[idx]/B[rowIdx--];

            if (new_t > t)
            {
                t = new_t;
                l = rowIdx;
            }
        }
    }

    //if ( unbounded )  // need to add check
}

int main (int argc, char *argv[])
{

    return 0;
}

