/***************************************************************/
/* This header file handles NRvector,NRmatrix,NRmat3d,Max, Min,*/
/*	SWAP and Random number generators for common use.		   */
/***************************************************************/

#ifndef _NR3_H_
#define _NR3_H_

//#define _CHECKBOUNDS_ 1
#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

#include <ctype.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

#define PREP(Fpre, Fnum, R_W)                    \
  {                                              \
    FILE *fp;                                    \
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

template <class T>
inline T SQR(const T a) {
  return a * a;
}

template <class T>
inline const T &MAX(const T &a, const T &b) {
  return b > a ? (b) : (a);
}

inline double MAX(const double &a, const double &b) {
  return (double)(b > a ? (b) : (a));
}

template <class T>
inline const T &MIN(const T &a, const T &b) {
  return b < a ? (b) : (a);
}

inline double MIN(const double &a, const double &b) {
  return (double)(b < a ? (b) : (a));
}

template <class T>
inline T SIGN(const T &a, const T &b) {
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline double SIGN(const double &a, const double &b) {
  return (double)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

template <class T>
inline void SWAP(T &a, T &b) {
  T dum = a;
  a = b;
  b = dum;
}

// exception handling

#ifndef _USENRERRORCLASS_
#define throw(message){printf("ERROR: %s\n     in file %s at line %d\n", \
                              message, __FILE__, __LINE__);              \
                       throw(1); }
#else
struct NRerror {
  char *message;
  char *file;
  int line;
  NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message, __FILE__, __LINE__));
void NRcatch(NRerror err) {
  printf("ERROR: %s\n     in file %s at line %d\n", err.message, err.file,
         err.line);
  exit(1);
}
#endif

// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)

// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define NRvector vector
#else

template <class T>
class NRvector {
 private:
  int nn;  // size of array. upper index is nn-1
  T *v;

 public:
  NRvector();
  explicit NRvector(int n);                  // Zero-based array
  NRvector(int n, const T &a);               // initialize to constant value
  NRvector(int n, const T *a);               // Initialize to array
  NRvector(const NRvector &rhs);             // Copy constructor
  NRvector &operator=(const NRvector &rhs);  // assignment
  typedef T value_type;                      // make T available externally
  T &operator[](const int i);                // i'th element
  const T &operator[](const int i) const;
  int size() const;
  void resize(int newn);              // resize (contents not preserved)
  void assign(int newn, const T &a);  // resize and assign a constant value
  ~NRvector();
};

/* NRvector definitions. */

template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), v(n > 0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(int n, const T &a) : nn(n), v(n > 0 ? new T[n] : NULL) {
  for (int i = 0; i < n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n > 0 ? new T[n] : NULL) {
  for (int i = 0; i < n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs)
    : nn(rhs.nn), v(nn > 0 ? new T[nn] : NULL) {
  for (int i = 0; i < nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T> &NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
  if (this != &rhs) {
    if (nn != rhs.nn) {
      if (v != NULL) delete[](v);
      nn = rhs.nn;
      v = nn > 0 ? new T[nn] : NULL;
    }
    for (int i = 0; i < nn; i++) {
      v[i] = rhs[i];
    }
  }
  return *this;
}

template <class T>
T &NRvector<T>::operator[](const int i)  // subscripting
{
#ifdef _CHECKBOUNDS_
  if (i < 0 || i >= nn) {
    throw("NRvector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
const T &NRvector<T>::operator[](const int i) const  // subscripting
{
#ifdef _CHECKBOUNDS_
  if (i < 0 || i >= nn) {
    throw("NRvector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
int NRvector<T>::size() const {
  return nn;
}

template <class T>
void NRvector<T>::resize(int newn) {
  if (newn != nn) {
    if (v != NULL) delete[](v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
}

template <class T>
void NRvector<T>::assign(int newn, const T &a) {
  if (newn != nn) {
    if (v != NULL) delete[](v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  for (int i = 0; i < nn; i++) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector() {
  if (v != NULL) delete[](v);
}

/* End of NRvector definitions. */

#endif  // ifdef _USESTDVECTOR_

template <class T>
class NRmatrix {
 private:
  int nn;
  int mm;
  T **v;

 public:
  NRmatrix();
  NRmatrix(int n, int m);                    // Zero-based array
  NRmatrix(int n, int m, const T &a);        // Initialize to constant
  NRmatrix(int n, int m, const T *a);        // Initialize to array
  NRmatrix(const NRmatrix &rhs);             // Copy constructor
  NRmatrix &operator=(const NRmatrix &rhs);  // assignment
  typedef T value_type;                      // make T available externally
  T *operator[](const int i);                // subscripting: pointer to row i
  const T *operator[](const int i) const;
  int nrows() const;
  int ncols() const;
  void resize(int newn, int newm);  // resize (contents not preserved)
  void assign(int newn, int newm,
              const T &a);  // resize and assign a constant value
  ~NRmatrix();
};

/* NRmatrix definitions. */

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m)
    : nn(n), mm(m), v(n > 0 ? new T *[n] : NULL) {
  int i, nel = m * n;
  if (v) v[0] = nel > 0 ? new T[nel] : NULL;
  for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a)
    : nn(n), mm(m), v(n > 0 ? new T *[n] : NULL) {
  int i, j, nel = m * n;
  if (v) v[0] = nel > 0 ? new T[nel] : NULL;
  for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a)
    : nn(n), mm(m), v(n > 0 ? new T *[n] : NULL) {
  int i, j, nel = m * n;
  if (v) v[0] = nel > 0 ? new T[nel] : NULL;
  for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs)
    : nn(rhs.nn), mm(rhs.mm), v(nn > 0 ? new T *[nn] : NULL) {
  int i, j, nel = mm * nn;
  if (v) v[0] = nel > 0 ? new T[nel] : NULL;
  for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
  for (i = 0; i < nn; i++)
    for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> &NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
  if (this != &rhs) {
    int i, j, nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[](v[0]);
        delete[](v);
      }
      nn = rhs.nn;
      mm = rhs.mm;
      v = nn > 0 ? new T *[nn] : NULL;
      nel = mm * nn;
      if (v) v[0] = nel > 0 ? new T[nel] : NULL;
      for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
    }
    for (i = 0; i < nn; i++)
      for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}

template <class T>
T *NRmatrix<T>::operator[](const int i)  // subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
  if (i < 0 || i >= nn) {
    throw("NRmatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
const T *NRmatrix<T>::operator[](const int i) const {
#ifdef _CHECKBOUNDS_
  if (i < 0 || i >= nn) {
    throw("NRmatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
int NRmatrix<T>::nrows() const {
  return nn;
}

template <class T>
int NRmatrix<T>::ncols() const {
  return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm) {
  int i, nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[](v[0]);
      delete[](v);
    }
    nn = newn;
    mm = newm;
    v = nn > 0 ? new T *[nn] : NULL;
    nel = mm * nn;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
  }
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T &a) {
  int i, j, nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[](v[0]);
      delete[](v);
    }
    nn = newn;
    mm = newm;
    v = nn > 0 ? new T *[nn] : NULL;
    nel = mm * nn;
    if (v) v[0] = nel > 0 ? new T[nel] : NULL;
    for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
  }
  for (i = 0; i < nn; i++)
    for (j = 0; j < mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix() {
  if (v != NULL) {
    delete[](v[0]);
    delete[](v);
  }
}

/* End of NRmatrix definitions. */

template <class T>
class NRMat3d {
 private:
  int nn;
  int mm;
  int kk;
  T ***v;

 public:
  NRMat3d();
  NRMat3d(int n, int m, int k);
  T **operator[](const int i);  // subscripting: pointer to row i
  const T *const *operator[](const int i) const;
  int dim1() const;
  int dim2() const;
  int dim3() const;
  ~NRMat3d();
};

/* NRmat3d definitions. */

template <class T>
NRMat3d<T>::NRMat3d() : nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T **[n]) {
  int i, j;
  v[0] = new T *[n * m];
  v[0][0] = new T[n * m * k];
  for (j = 1; j < m; j++) v[0][j] = v[0][j - 1] + k;
  for (i = 1; i < n; i++) {
    v[i] = v[i - 1] + m;
    v[i][0] = v[i - 1][0] + m * k;
    for (j = 1; j < m; j++) v[i][j] = v[i][j - 1] + k;
  }
}

template <class T>
T **NRMat3d<T>::operator[](const int i)  // subscripting: pointer to row i
{
  return v[i];
}

template <class T>
const T *const *NRMat3d<T>::operator[](const int i) const {
  return v[i];
}

template <class T>
int NRMat3d<T>::dim1() const {
  return nn;
}

template <class T>
int NRMat3d<T>::dim2() const {
  return mm;
}

template <class T>
int NRMat3d<T>::dim3() const {
  return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d() {
  if (v != NULL) {
    delete[](v[0][0]);
    delete[](v[0]);
    delete[](v);
  }
}
/* End of NRmat3d definitions. */

/////////////////////////////////////////////////////////////////////////
// basic type names (redefine if your bit lengths don't match)

typedef int Int;  // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong;  // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong;  // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char;  // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub;  // default floating type
typedef long double Ldoub;

typedef complex<double> Complex;  // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

// Uint proto_nan[2]={0xffffffff, 0x7fffffff};
// double NaN = *( double* )proto_nan;

// Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char *> VecCharp_I;
typedef NRvector<Char *> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub *> VecDoubp_I;
typedef NRvector<Doub *> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

////////////////////////////////////////////////////////////////////
/*
        Classes for Random Number Generators Below.
*/

class Ran {
  /*
  Implementation of the highest quality recommended generator. The constructor
  is called with an integer seed and creates an instance of the generator. The
  member functions int64, doub, and int32 return the next values in the random
  sequence, as a variable type indicated by their names. The period of the
  generator is =3.138*10^57.
  */

 public:
  Ran(Ullong j) : v(4101842887655102017LL), w(1) {
    u = j ^ v;
    int64();
    v = u;
    int64();
    w = v;
    int64();
  }
  //	Constructor. Call with any integer seed (except value of v above).  Ran
  // myran(seed)

  ~Ran(){};

  /*Return 64-bit random integer. See text for explanation of method.*/
  inline Ullong int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);

    Ullong x = u ^ (u << 21);
    x ^= x >> 35;
    x ^= x << 4;

    return (x + v) ^ w;
  }

  /*Return random double-precision floating value in the range 0. to 1.*/
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  /*Return 32-bit random integer.*/
  inline Uint int32() { return (Uint)int64(); }

 private:
  Ullong u, v, w;
};

class Ranq1 {
  /*
  Recommended generator for everyday use. The period is 1.8*10^19. Calling
  conventions same as Ran, above.
  */

 public:
  Ranq1(Ullong j) : v(4101842887655102017LL) {
    v ^= j;
    v = int64();
  }

  ~Ranq1(){};

  inline Ullong int64() {
    v ^= v >> 21;
    v ^= v << 35;
    v ^= v >> 4;
    return v * 2685821657736338717LL;
  }

  inline double doub() { return 5.42101086242752217E-20 * int64(); }

  inline Uint int32() { return (Uint)int64(); }

 private:
  Ullong v;
};

class Ranq2 {
  /*
  Backup generator if Ranq1 has too short a period and Ran is too slow. The
  period is 8.5*10^37. Calling conventions same as Ran, above.
  */

 public:
  Ranq2(Ullong j) : v(4101842887655102017LL), w(1) {
    v ^= j;
    w = int64();
    v = int64();
  }
  ~Ranq2(){};

  inline Ullong int64() {
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    return v ^ w;
  }

  inline double doub() { return 5.42101086242752217E-20 * int64(); }

  inline Uint int32() { return (Uint)int64(); }

 private:
  Ullong v, w;
};

class Ranfib {
  /*
  Implements Knuth's subtractive generator using only floating operations. See
  text for cautions.
  */

 public:
  //	Constructor. Call with any integer seed. Uses Ranq1 to initialize.
  Ranfib(Ullong j) : inext(0), inextp(31) {
    Ranq1 init(j);
    for (int k = 0; k < 55; k++) {
      dtab[k] = init.doub();
    }
  }
  ~Ranfib(){};

  // Returns random double-precision floating value between 0. and 1.
  double doub() {
    if (++inext == 55) inext = 0;
    if (++inextp == 55) inextp = 0;
    dd = dtab[inext] - dtab[inextp];
    if (dd < 0) dd += 1.0;
    return (dtab[inext] = dd);
  }

  //	Returns random 32-bit integer. Recommended only for testing purposes.
  inline unsigned long int32() {
    return (unsigned long)(doub() * 4294967295.0);
  }

 private:
  double dtab[55];
  double dd;

  int inext;
  int inextp;
};

/*Normal deviates' class.*/
class Normaldev_BM : public Ran {
 public:
  // Constructor arguments are mu,sig and a random sequence seed.
  Normaldev_BM(double mmu, double ssig, Ullong i)
      : Ran(i), mu(mmu), sig(ssig), storedval(0.){};
  ~Normaldev_BM(){};

  // Return a normal deviate.
  double dev() {
    double v1, v2, rsq, fac;
    if (storedval == 0.) {  // We don't have an extra deviate handy, so
      do {
        v1 =
            2.0 * doub() - 1.0;  // pick two uniform numbers in the square exv2=
        v2 = 2.0 * doub() - 1.0;  // tending from -1 to +1 in each direction,
        rsq = v1 * v1 +
              v2 * v2;  // see if they are in the unit circle, or try again.
      } while (rsq >= 1.0 || rsq == 0.0);

      fac = sqrt(-2.0 * log(rsq) /
                 rsq);  // Now make the Box-Muller transformation to
                        // get two normal deviates. Return one and
                        // save the other for next time.
      storedval = v1 * fac;
      return mu + sig * v2 * fac;
    } else {  // We have an extra deviate handy,
      fac = storedval;
      storedval = 0.;
      return mu + sig * fac;  // so return it.
    }
  }

 private:
  double mu;
  double sig;
  double storedval;
};

/*Normal deviates.*/
class Normaldev : public Ran {
 public:
  Normaldev(double mmu, double ssig, Ullong i) : Ran(i), mu(mmu), sig(ssig){};
  ~Normaldev(){};

  double dev() {
    double u, v, x, y, q;
    do {
      u = doub();
      v = 1.7156 * (doub() - 0.5);
      x = u - 0.449871;
      y = abs(v) + 0.386595;
      q = SQR(x) + y * (0.19600 * y - 0.25472 * x);
    } while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4. * log(u) * SQR(u)));

    return mu + sig * v / u;
  }

 private:
  double mu;
  double sig;
};

//! Solves for a vector u[0..n-1] the tridiagonal linear set given by equation
//! (2.4.1). a[0..n-1],
// b[0..n-1], c[0..n-1], and r[0..n-1] are input vectors and are not modified.
void tridag(VecDoub_I &a, VecDoub_I &b, VecDoub_I &c, VecDoub_I &r,
            VecDoub_O &u) {
  int j, n = static_cast<int>(a.size());
  double bet;
  VecDoub gam(n);  // One vector of workspace, gam, is needed.
  if (b[0] == 0.0)
    throw("Error 1 in tridag");  // If this happens, then you should rewrite
                                 // your equations as a set of order N - 1, with
                                 // u1 trivially eliminated.
  u[0] = r[0] / (bet = b[0]);
  for (j = 1; j < n; j++)  // Decomposition and forward substitution.
  {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) throw("Error 2 in tridag");  // Algorithm fails; see below.
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (j = n - 2; j >= 0; j--)
    u[j] -= gam[j + 1] * u[j + 1];  // Backsubstitution.
}

//! Solves for a vector x[0..n-1] the ¡°cyclic¡± set of linear equations given by
//! equation (2.7.9).
// a,b, c, and r are input vectors, all dimensioned as [0..n-1],
// while alpha and beta are the corner entries in the matrix.
// The input is not modified.
void cyclic(VecDoub_I &a, VecDoub_I &b, VecDoub_I &c, const Doub alpha,
            const Doub beta, VecDoub_I &r, VecDoub_O &x) {
  int i, n = static_cast<int>(a.size());
  double fact, gamma;
  if (n <= 2) throw("n too small in cyclic");
  VecDoub bb(n), u(n), z(n);
  gamma = -b[0];  // Avoid subtraction error in forming bb[0].
  bb[0] =
      b[0] - gamma;  // Set up the diagonal of the modified tridiagonal system.
  bb[n - 1] = bb[n - 1] - alpha * beta / gamma;
  for (i = 1; i < n - 1; i++) bb[i] = b[i];
  tridag(a, bb, c, r, x);  // Solve A.x =r.
  u[0] = gamma;            // Set up the vector u.
  u[n - 1] = alpha;
  for (i = 1; i < n - 1; i++) u[i] = 0.0;
  tridag(a, bb, c, u, z);  // Solve A.z = u.
  fact = (x[0] + beta * x[n - 1] / gamma) /
         (1.0 + z[0] + beta * z[n - 1] / gamma);  // Form v.x/(1+v.z).
  for (i = 0; i < n; i++)
    x[i] -= fact * z[i];  // Now get the solution vector x.
}

// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
  turn_on_floating_exceptions() {
    int cw = _controlfp(0, 0);
    cw &= ~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE);
    _controlfp(cw, MCW_EM);
  }
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

#endif /* _NR3_H_ */
