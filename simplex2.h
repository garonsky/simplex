int getIndex(int nRows, int nCols, int r, int c)
{
    r--;
    c--;

    assert(r>=0 && r < nRows);
    assert(c>=0 && c < nCols);
    return c*nRows+r;
}


void printArray(int m, int n, double a[], int transpose)
{
    int rowC = 0;
    int colC = 0;

    while (rowC < m)
    {
        while(colC < n)
        {
            if (transpose == 0)
            {
                printf (" a[%d,%d] = %f\t", rowC+1, colC+1, a[getIndex(m, n, rowC+1, colC+1)]);
            }
            else
            {
                printf (" a[%d,%d] = %f\n", rowC+1, colC+1, a[getIndex(m, n, rowC+1, colC+1)]);
            }
            colC++;
        }
        printf("\n");
        colC=0;
        rowC++;
    }
}

void printMatrix(double *X, int m, int n);

void printMatrix2(double *An, int m, int n, double *XB, double *Xb, double *c, double *cb, double t, int e, int l);

void printLine(int nFactor);

void copyRow(double *d, const double* An, int m, int n, int l);
void copyColumn(double *r, const double* An, int m, int n, int e, double a);
void addToColumn(double *An, int m, int n, int e, const double *t);
void subtractFromColumn(double *An, int m, int n, int e, const double* t, double *XB, double *Xb);

// c has a length of n
//int argmax(const double *An, int m, int n, const double *c);

double calcEuclideanDistance(const double *d, int m, int n, int c);

//int findEnteringVariable(const double *An, int m, int n, const double *c);

int findLeavingVariable(const double *b, const double *An, int m, int n, int e, int *l, double *t);

int pivot(const double *An, int m, int n);

//void expand(const double *B, const double *An, int m, int n, int e, int *l, double *t);

void swap(double *An, int m, int n, int c1, double *XB, int c2);

