int getIndex(int nRows, int nCols, int r, int c);
void printArray(int m, int n, double *c, int transpose, const char* label);
void printMatrix(double *a, int m, int n, double *b, double *c, double z, const char *label, int pl, int pe);
void printLine(int nFactor);
double calcEuclideanDistance(const double *d, int m, int n, int c);

int findEnteringVariable(const double *An, int m, int n, const double *c);
int findLeavingVariable(const double *b, const double *An, int m, int n, int e, int *l, double *t);
void eliminate(double *a, int m, int n, int e, int l, double *b, double *c, double *z);
void copyRow(double *r, const double* An, int m, int n, int l);
void copyRow2(const double *r, double* An, int m, int n, int l);

void swapColumn(double *An, int m, int n, int c1, double *XB, int c2);


#ifdef false 

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
   a[m*5] = 0;
   a[m*5+1] = 0;
   a[m*5+2] = 0;
   a[m*5+2] = 1;

#endif


#ifdef false
   b[0] = 29;
   b[1] = -10;
   b[2] = 3;
   b[3] = 20;
   b[4] = 20;

   c[0] = -2;
   c[1] = -4;
   c[2] = 4;
   c[3] = -9;
   c[4] = 0;
   c[5] = 0;
   c[6] = 0;
   c[7] = 0;
   c[8] = 0;

   /* The elements of the first column */
   a[0] = 2;
   a[1] = -1;
   a[2] = 0;
   a[3] = 1;
   a[4] = 0;
   /* The elements of the second column */
   a[m] = 3;
   a[m+1] = -1;
   a[m+2] = 1;
   a[m+3] = 1;
   a[m+4] = 0;
   /* The elements of the third column */
   a[m*2] = 0;
   a[m*2+1] = 0;
   a[m*2+2] = 0;
   a[m*2+3] = 0;
   a[m*2+4] = 1;
   /* The elements of the fourth column */
   a[m*3]   = 0;
   a[m*3+1] = 0;
   a[m*3+2] = 0;
   a[m*3+3] = -1;
   a[m*3+4] = 1;

   /* The elements of the fifth column */
   a[m*4]   = 1;
   a[m*4+1] = 0;
   a[m*4+2] = 0;
   a[m*4+3] = 0;
   a[m*4+4] = 0;


   /* The elements of the sixth column */
   a[m*5]   = 0;
   a[m*5+1] = 1;
   a[m*5+2] = 0;
   a[m*5+3] = 0;
   a[m*5+4] = 0;

   /* The elements of the seventh column */
   a[m*6]   = 0;
   a[m*6+1] = 0;
   a[m*6+2] = 1;
   a[m*6+3] = 0;
   a[m*6+4] = 0;


   /* The elements of the eighth column */
   a[m*7]   = 0;
   a[m*7+1] = 0;
   a[m*7+2] = 0;
   a[m*7+3] = 1;
   a[m*7+4] = 0;


   /* The elements of the ninth column */
   a[m*8]   = 0;
   a[m*8+1] = 0;
   a[m*8+2] = 0;
   a[m*8+3] = 0;
   a[m*8+4] = 1;

   //x1=10; x2=3.0; x3=0; x4=20; z = 212

#endif

#ifdef false

// TEST SAMPLE

  m = 3; /* Size of Column ( the number of rows ) */
  n = 6; /* Size of Row ( the number of columns ) */
  lda = 3; /* Leading dimension of 5 * 4 matrix is 5 */

  incx = 1;
  incy = 1;
  alpha = 1;
  beta = 0;

  a = (double *)malloc(sizeof(double)*m*(n+m));
  b = (double *)malloc(sizeof(double)*n);
  c = (double *)malloc(sizeof(double)*(n+m));

  b[0] = 14;
  b[1] = 28;
  b[2] = 30;

  c[0] = -1;
  c[1] = -2;
  c[2] = 1;
  c[3] = 0;
  c[4] = 0;
  c[5] = 0;

  /* The elements of the first column */
  a[0] = 2;
  a[1] = 4;
  a[2] = 2;

  /* The elements of the second column */
  a[m] = 1;
  a[m+1] = 2;
  a[m+2] = 5;

  /* The elements of the third column */
  a[m*2] = 1;
  a[m*2+1] = 3;
  a[m*2+2] = 5;

  /* The elements of the fourth column */
  a[m*3] = 1;
  a[m*3+1] = 0;
  a[m*3+2] = 0;


  /* The elements of the fifth column */
  a[m*4] = 0;
  a[m*4+1] = 1;
  a[m*4+2] = 1;


  /* The elements of the sixth column */
  a[m*5] = 0;
  a[m*5+1] = 0;
  a[m*5+2] = 1;

  //x1=5; x2=4; z = 13

#endif


#ifdef false

// TEST SAMPLE

  m = 3; /* Size of Column ( the number of rows ) */
  n = 7; /* Size of Row ( the number of columns ) */
  lda = 4; /* Leading dimension of 5 * 4 matrix is 5 */

  incx = 1;
  incy = 1;
  alpha = 1;
  beta = 0;

  a = (double *)malloc(sizeof(double)*m*(n+m));
  b = (double *)malloc(sizeof(double)*n);
  c = (double *)malloc(sizeof(double)*(n+m));

  b[0] = 12;
  b[1] = 7;
  b[2] = 10;

  c[0] = -2;
  c[1] = -4;
  c[2] = -3;
  c[3] = -1;
  c[4] = 0;
  c[5] = 0;
  c[6] = 0;

  /* The elements of the first column */
  a[0] = 3;
  a[1] = 1;
  a[2] = 2;

  /* The elements of the second column */
  a[m]   = 1;
  a[m+1] = -3;
  a[m+2] = 1;

  /* The elements of the third column */
  a[m*2]   = 1;
  a[m*2+1] = 2;
  a[m*2+2] = 3;

  /* The elements of the fourth column */
  a[m*3]   = 4;
  a[m*3+1] = 3;
  a[m*3+2] = -1;


  /* The elements of the fifth column */
  a[m*4]   = 1;
  a[m*4+1] = 0;
  a[m*4+2] = 0;


  /* The elements of the sixth column */
  a[m*5]   = 0;
  a[m*5+1] = 1;
  a[m*5+2] = 0;

  /* The elements of the sixth column */
  a[m*6]   = 0;
  a[m*6+1] = 0;
  a[m*6+2] = 1;

  //optimal: x1 = 0; x2 = 10.4; x3 = 0; x4 = 0.4; x5 = 0; x6 = 37.0; x7 = 0;

#endif
