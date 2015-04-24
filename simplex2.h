int getIndex(int nRows, int nCols, int r, int c)
{
    r--;
    c--;

    assert(r>=0 && r < nRows);
    assert(c>=0 && c < nCols);
    return c*nRows+r;
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

void printMatrix(double *X, int m, int n);

void copyRow(double *d, const double* An, int m, int n, int l);
void copyColumn(double *r, const double* An, int m, int n, int e);
void addToColumn(const double *An, int m, int n, int e, const double *t);
void addToColumn2(const double *An, int m, int n, int e, const double* t, double *XB, int *Xb);

// c has a length of n
//int argmax(const double *An, int m, int n, const double *c);

double calcEuclideanDistance(const double *d, int m, int n, int c);

//int findEnteringVariable(const double *An, int m, int n, const double *c);

int findLeavingVariable(const double *b, const double *An, int m, int n, int e, int *l, double *t);

int pivot(const double *An, int m, int n);

//void expand(const double *B, const double *An, int m, int n, int e, int *l, double *t);

