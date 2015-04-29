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
