int getIndex(int nRows, int nCols, int r, int c);
void printArray(int m, int n, double *c, int transpose, const char* label);
void printMatrix(double *a, int m, int n, double *b, double *c, double z, const char *label);
void printLine(int nFactor);
double calcEuclideanDistance(const double *d, int m, int n, int c);

int findEnteringVariable(const double *An, int m, int n, const double *c);
int findEnteringVariable2(const double *An, int m, int n, const double *c);
int findEnteringVariable3(const double *An, int m, int n, const double *c);
int findLeavingVariable(const double *b, const double *An, int m, int n, int e, int *l, double *t);
void eliminate(double *a, int m, int n, int e, int l, double *b, double *c, double *z);
void copyRow(double *r, const double* An, int m, int n, int l);
void copyRow2(const double *r, double* An, int m, int n, int l);

void swapColumn(double *An, int m, int n, int c1, double *XB, int c2);


#ifdef 0

// TEST SAMPLE

   m = 4; /* Size of Column ( the number of rows ) */
   n = 6; /* Size of Row ( the number of columns ) */
   lda = 4; /* Leading dimension of 5 * 4 matrix is 5 */

   incx = 1;
   incy = 1;
   alpha = 1;
   beta = 0;

   a = (double *)malloc(sizeof(double)*m*(n+m));
   b = (double *)malloc(sizeof(double)*n);
   c = (double *)malloc(sizeof(double)*(n+m));

   b[0] = 6;
   b[1] = 3;
   b[2] = 5;
   b[3] = 4;

   c[0] = -4;
   c[1] = -3;
   c[2] = 0;
   c[3] = 0;
   c[4] = 0;
   c[5] = 0;

   /* The elements of the first column */
   a[0] = 2;
   a[1] = -3;
   a[2] = 0;
   a[3] = 2;
   /* The elements of the second column */
   a[m] = 3;
   a[m+1] = 2;
   a[m+2] = 2;
   a[m+3] = 1;
   /* The elements of the third column */
   a[m*2] = 1;
   a[m*2+1] = 0;
   a[m*2+2] = 0;
   a[m*2+3] = 0;
   /* The elements of the fourth column */
   a[m*3] = 0;
   a[m*3+1] = 1;
   a[m*3+2] = 0;
   a[m*3+3] = 0;

   /* The elements of the fifth column */
   a[m*4] = 0;
   a[m*4+1] = 0;
   a[m*4+2] = 1;
   a[m*4+2] = 0;

   /* The elements of the sixth column */
   a[m*4] = 0;
   a[m*4+1] = 0;
   a[m*4+2] = 0;
   a[m*4+2] = 1;

#endif
