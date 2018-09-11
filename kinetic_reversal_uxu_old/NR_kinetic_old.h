#ifndef NR_KINETIC_OLD_H
#define NR_KINETIC_OLD_H

#include <stddef.h>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <ctime>

#define PREP(Fpre, Fnum, R_W)                    \
  {                                              \
    FILE* fp;                                    \
    errno_t err;                                 \
    {                                            \
      char PRT[100];                             \
      sprintf_s(PRT, "%s%dsA.dat", Fpre, Fnum);  \
      if ((err = fopen_s(&fp, PRT, R_W)) != 0) { \
        printf_s("Fale File=%s.\n", PRT);        \
        exit(0);                                 \
      }                                          \
    }
#define CLSP      \
  { fclose(fp); } \
  }

#define NR_END 1
#define FREE_ARG char*

template <class T>
inline T SQR(const T a) {
  return a * a;
}

void nrerror(const char error_text[])
{
  fprintf(stderr, "Numerial Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

int* ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int* v;
  v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));

  if (!v) nrerror("allocation failure in ivector()");

  return v - nl + NR_END;
}

void free_ivector(int* v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG)(v + nl - NR_END));
}

double* vector(long nl, long nh) {
  double* v;
  v = (double*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));

  if (!v) nrerror("allocation failure in vector()");

  return v - nl + NR_END;
}

void free_vector(double* v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
}

struct bdst** adst_vector(long nl, long nh) {
  struct bdst** v;
  v = (struct bdst**)malloc(
      (size_t)((nh - nl + 1 + NR_END) * sizeof(struct bdst*)));

  if (!v) nrerror("allocation failure in adst_vector()");

  return v - nl + NR_END;
}

struct bdln** adln_vector(long nl, long nh) {
  struct bdln** v;
  v = (struct bdln**)malloc(
      (size_t)((nh - nl + 1 + NR_END) * sizeof(struct bdln*)));

  if (!v) nrerror("allocation failure in adst_vector()");

  return v - nl + NR_END;
}

void free_adst_vector(struct bdst** v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
}

void free_adln_vector(struct bdln** v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
}

int* int_vector(long nl, long nh) {
  int* v;
  v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));

  if (!v) nrerror("allocation failure in vector()");

  return v - nl + NR_END;
}

void free_int_vector(int* v, long nl, long nh) {
  free((FREE_ARG)(v + nl - NR_END));
}

double** matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double** m;
  m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));

  if (!m) nrerror("allocation failure 1 in matrix()");

  m += NR_END;
  m -= nrl;
  m[nrl] = (double*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));

  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  return m;
}

float** float_matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float** m;
  m = (float**)malloc((size_t)((nrow + NR_END) * sizeof(float*)));

  if (!m) nrerror("allocation failure 1 in matrix()");

  m += NR_END;
  m -= nrl;
  m[nrl] = (float*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));

  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  return m;
}

int** int_matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int** m;
  m = (int**)malloc((size_t)((nrow + NR_END) * sizeof(int*)));

  if (!m) nrerror("allocation failure 1 in matrix()");

  m += NR_END;
  m -= nrl;
  m[nrl] = (int*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));

  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");

  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

  return m;
}

void free_int_matrix(int** m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_matrix(double** m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void free_float_matrix(float** m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

double*** f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) {
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  double*** t;
  t = (double***)malloc((size_t)((nrow + NR_END) * sizeof(double**)));

  if (!t) nrerror("allocation failure 1 in f3tensor()");

  t += NR_END;
  t -= nrl;
  t[nrl] = (double**)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double*)));

  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");

  t[nrl] += NR_END;
  t[nrl] -= ncl;
  t[nrl][ncl] =
      (double*)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(double)));

  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");

  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;

  for (i = nrl + 1; i <= nrh; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;

    for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
  }

  return t;
}

void free_f3tensor(double*** t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh) {
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
}

void tridag(double a[], double b[], double c[], double r[], double u[],
            unsigned long n)
/*Solves for a vector u[1..n] the tridiagonal linear set given by equation
(2.4.1). a[1..n], b[1..n], c[1..n], and r[1..n] are input vectors and are not
modied.*/
{
  unsigned long j;
  double bet, *gam;
  gam = vector(1, n); /* One vector of workspace, gam is needed.*/

  if (b[1] == 0.0) nrerror("Error 1 in tridag");

  /*If this happens then you should rewrite your equations as a set of order
  N-1, with u2 trivially eliminated.*/
  u[1] = r[1] / (bet = b[1]);

  for (j = 2; j <= n; j++) /*Decomposition and forward substitution.*/
  {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];

    if (bet == 0.0)
      nrerror("Error 2 in tridag"); /* Algorithm fails; see below*/

    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }

  for (j = (n - 1); j >= 1; j--)
    u[j] -= gam[j + 1] * u[j + 1]; /*Backsubstitution.*/

  free_vector(gam, 1, n);
}

void cyclic(double a[], double b[], double c[], double alpha, double beta,
            double r[], double x[], unsigned long n)
/*Solves for a vector x[1..n] the \cyclic" set of linear equations given by
equation (2.7.9). a, b, c, and r are input vectors, all dimensioned as [1..n],
while alpha and beta are the corner
entries in the matrix. The input is not modied.*/
{
  void tridag(double a[], double b[], double c[], double r[], double u[],
              unsigned long n);
  unsigned long i;
  double fact, gamma, *bb, *u, *z;

  if (n <= 2) nrerror("n too small in cyclic");

  bb = vector(1, n);
  u = vector(1, n);
  z = vector(1, n);
  gamma = -b[1];        /* Avoid subtraction error in forming bb[1].*/
  bb[1] = b[1] - gamma; /*Set up the diagonal of the modied tridiagonalsystem*/
  bb[n] = b[n] - alpha * beta / gamma;

  for (i = 2; i < n; i++) bb[i] = b[i];

  tridag(a, bb, c, r, x, n); /* Solve A . x = r.*/
  u[1] = gamma;              /*Set up the vector u.*/
  u[n] = alpha;

  for (i = 2; i < n; i++) u[i] = 0.0;

  tridag(a, bb, c, u, z, n); /*Solve A . z = u.*/
  fact = (x[1] + beta * x[n] / gamma) /
         (1.0 + z[1] + beta * z[n] / gamma); /*Form v . x=(1 + v . z).*/

  for (i = 1; i <= n; i++)
    x[i] -= fact * z[i]; /*Now get the solution vector x.*/

  free_vector(z, 1, n);
  free_vector(u, 1, n);
  free_vector(bb, 1, n);
}

#define IAR1 16807
#define IMR1 2147483647
#define AMR1 (1.0 / IMR1)
#define IQR1 127773
#define IRR1 2836
#define NTABR1 32
#define NDIVR1 (1 + (IMR1 - 1) / NTABR1)
#define EPSR1 1.2e-7
#define RNMXR1 (1.0 - EPSR1)

double ran2(long* idum)
/*\Minimal" random number generator of Park and Miller with Bays-Durham shu.e
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0
(exclusive of the endpoint values). Call with idum a negative integer to
initialize; thereafter, do not alter idum between successive deviates in a
sequence. RNMX should approximate the largest doubleing value that is less
than 1.*/
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTABR1];
  double temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1)
      *idum = 1;
    else
      *idum = -(*idum);

    for (j = NTABR1 + 7; j >= 0; j--) {
      k = (*idum) / IQR1;
      *idum = IAR1 * (*idum - k * IQR1) - IRR1 * k;

      if (*idum < 0) *idum += IMR1;

      if (j < NTABR1) iv[j] = *idum;
    }

    iy = iv[0];
  }

  k = (*idum) / IQR1;
  *idum = IAR1 * (*idum - k * IQR1) - IRR1 * k;

  if (*idum < 0) *idum += IMR1;

  j = iy / NDIVR1;
  iy = iv[j];
  iv[j] = *idum;

  if ((temp = AMR1 * iy) > RNMXR1)
    return RNMXR1;
  else
    return temp;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0 / MBIG)
/*According to Knuth, any large MBIG, and any smaller (but still large) MSEED
can be substituted for the above values.*/
double ran1(long* idum)
/*Returns a uniform random deviate between 0:0 and 1:0. Set idum to any negative
value to initialize or reinitialize the sequence.*/
{
  static int inext, inextp;
  static long ma[56];
  static int iff = 0;
  /*The value 56 (range ma[1..55]) is special and
  should not be modied; see Knuth.
  */
  long mj, mk;
  int i, ii, k;

  if (*idum < 0 || iff == 0)  // Initialization.
  {
    iff = 1;
    mj = MSEED -
         (*idum < 0 ? -*idum : *idum);  // Initialize ma[55] using the seed idum
                                        // and the large numberMSEED.
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;

    for (i = 1; i <= 54; i++)  // Now initialize the rest of the table,
    {
      ii = (21 * i) % 55;  // in a slightly random order,
      ma[ii] = mk;         // with numbers that are not especially random.
      mk = mj - mk;

      if (mk < MZ) mk += MBIG;

      mj = ma[ii];
    }

    for (k = 1; k <= 4;
         k++)  // We randomize them by \warming up the generator."
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];

        if (ma[i] < MZ) ma[i] += MBIG;
      }

    inext = 0;    // Prepare indices for our rst generated number.
    inextp = 31;  // The constant 31 is special; see Knuth.
    *idum = 1;
  }

  // Here is where we start, except on initialization.
  if (++inext == 56) inext = 1;  // Increment inext and inextp, wrapping around

  if (++inextp == 56) inextp = 1;  // 56 to 1.

  mj = ma[inext] - ma[inextp];  // Generate a new random number subtractively.

  if (mj < MZ) mj += MBIG;  // Be sure that it is in range.

  ma[inext] = mj;   // Store it,
  return mj * FAC;  // and output the derived uniform deviate.
}
#endif  // !NR_KINETIC_OLD_H
