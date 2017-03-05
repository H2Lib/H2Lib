
/* ------------------------------------------------------------
 * This is the file "basic.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file basic.h
 *  @author Steffen B&ouml;rm */

#ifndef BASIC_H
#define BASIC_H

/** @defgroup basic basic
 *  @brief Miscellaneous auxiliary functions and macros.
 *  @{ */

/** @brief Stop watch for run-time measurements. */
typedef struct _stopwatch stopwatch;

/** @brief Pointer to a @ref stopwatch object. */
typedef stopwatch *pstopwatch;

#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "settings.h"

#ifdef USE_SIMD
#include "simd.h"
#endif

/* ------------------------------------------------------------
 * Miscellaneous
 * ------------------------------------------------------------ */

/** @brief Mathematical constant @f$\pi@f$. */
#ifndef M_PI
#define M_PI 3.141592653589793238462643383
#endif

/** @brief Number of power iteration steps used to approximate the spectral norm. */
#define NORM_STEPS 20

/** @brief Deflation level for NetCDF compression.
 *
 *  The value 9 corresponds to maximal compression. */
#define NC_DEFLATE_LEVEL 9

/** @brief Reasonable cut-off depth for parallelization. */
extern int max_pardepth;

/** @brief "Machine accuracy" for some algorithms */
#define H2_MACH_EPS 1e-13

/* ------------------------------------------------------------
 Set up the library
 ------------------------------------------------------------ */

/** @brief Initialize the library.
 *
 *  This function prepares the run-time environment for calls to
 *  library functions, e.g., by initializing external libraries
 *  like GTK+ or FreeGLUT that are required by some functions.
 *
 *  @param argc Number of command line parameters.
 *  @param argv Values of command line parameters.
 */
HEADER_PREFIX void
init_h2lib(int *argc, char ***argv);

/** @brief Uninitialize the library.
 *
 *  This function cleans up the run-time environment once the
 *  library is no longer neede. */
HEADER_PREFIX void
uninit_h2lib();

/* ------------------------------------------------------------
 * General utility macros and functions
 * ------------------------------------------------------------ */

/**
 * @brief Wrapper for either aligned or not align memory allocation function
 * @param sz Number of bytes that should be allocated.
 */
#ifdef USE_SIMD
#define h2_malloc(sz) _mm_malloc(sz, VALIGN)
#else
#define h2_malloc(sz) malloc(sz)
#endif

/**
 * @brief Wrapper for either aligned or not align memory free function
 * @param p Pointer to an object that should be freed.
 */
#ifdef USE_SIMD
#define h2_free(p) _mm_free(p)
#else
#define h2_free(p) free(p)
#endif


/**
 * @brief Current amount of allocated memory (bytes)
 * @return Current amount of currently allocated memory is returned.
 */
HEADER_PREFIX size_t
get_current_memory();

/**
 * @brief Round up some value x to the nearest multiple of N.
 * @param x Value that should be round up to nearest multiple of N.
 * @param N Muliples of N are used for rounding up x.
 * @return ((x + N - 1) / N) * N
 */
#define ROUNDUP(x, N) ((((x) + (N) - 1) / (N)) * (N))

#ifdef ABS
#undef ABS
#endif

#ifdef UINT_MAX
#undef UINT_MAX
#endif

/* Macros for the type "field" */

/** @brief Get the real part @f$a@f$ of a field element @f$x=a+ib@f$. */
#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define REAL(x) crealf(x)
#else
#define REAL(x) creal(x)
#endif
#else
#define REAL(x) (x)
#endif

/** @brief Get the imaginary part @f$b@f$ of a field element @f$x=a+ib@f$. */
#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define IMAG(x) cimagf(x)
#else
#define IMAG(x) cimag(x)
#endif
#else
#define IMAG(x) 0.0
#endif

/** @brief Compute the complex conjugate @f$\bar x@f$ of a field element @f$x@f$. */
#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define CONJ(x) conjf(x)
#else
#define CONJ(x) conj(x)
#endif
#else
#define CONJ(x) (x)
#endif

/** @brief Compute the square of the absolute value @f$|x|^2@f$ of a field element @f$x@f$. */
#define ABSSQR(x) _h2_abssqr(x)

/** @brief Compute the absolute value @f$|x|@f$ of a field element @f$x@f$. */
#define ABS(x) REAL_SQRT(ABSSQR(x))

/** @brief Compute the sign @f$\mathop{\rm sgn}(x)@f$ of a field
 *  element @f$x@f$ with the convention @f$\mathop{\rm sgn}(0)=1@f$. */
#define SIGN1(x) _h2_sgn1(x)

/** @brief Compute a (pseudo-)random field element */
#define FIELD_RAND() _h2_fieldrand()

/* Macros for the type "real" */

/** @brief Compute the absolute value @f$|x|@f$ of a real number @f$x@f$. */
#define REAL_ABS(x) fabs(x)

/** @brief Compute the square @f$x^2@f$ of a real number @f$x@f$. */
#define REAL_SQR(x) _h2_real_sqr(x)

/** @brief Compute the sine @f$\sin(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_SIN(x) sinf(x)
#else
#define REAL_SIN(x) sin(x)
#endif

/** @brief Compute the cosine @f$\cos(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_COS(x) cosf(x)
#else
#define REAL_COS(x) cos(x)
#endif

/** @brief Compute the arcsine @f$\asin(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_ASIN(x) asinf(x)
#else
#define REAL_ASIN(x) asin(x)
#endif

/** @brief Compute the arccosine @f$\acos(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_ACOS(x) acosf(x)
#else
#define REAL_ACOS(x) acos(x)
#endif

/** @brief Compute the tangent @f$\tan(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_TAN(x) tanf(x)
#else
#define REAL_TAN(x) tan(x)
#endif

/** @brief Compute the square root @f$\sqrt{x}@f$ of a non-negative real nunber @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_SQRT(x) sqrtf(x)
#else
#define REAL_SQRT(x) sqrt(x)
#endif

/** @brief Compute the reciprocal square root @f$1.0/\sqrt{x}@f$ of a
 *  non-negative real nunber @f$x@f$. */
#define REAL_RSQRT(x) _h2_rsqrt(x)

/** @brief Compute the @f$y@f$-th power @f$x^y@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_POW(x, y) powf(x, y)
#else
#define REAL_POW(x, y) pow(x, y)
#endif

/** @brief Compute the natural logarithm @f$\ln(x)@f$ of a positive real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_LOG(x) logf(x)
#else
#define REAL_LOG(x) log(x)
#endif

/** @brief Compute the exponetial function @f$\exp(x)@f$ of a real number @f$x@f$. */
#ifdef USE_FLOAT
#define REAL_EXP(x) expf(x)
#else
#define REAL_EXP(x) exp(x)
#endif

/** @brief Compute a (pseudo-)random real number */
#define REAL_RAND() _h2_realrand()

/** @brief Compute dot product of vectors of dimension 2, i.e.
 *   @f$ \bar x_1 y_1 + \bar x_2 y_2 @f$ */
#define DOT2(x,y) (CONJ((x)[0]) * (y)[0] + CONJ((x)[1]) * (y)[1])

/** @brief Compute dot product of vectors of dimension 3, i.e.
 *   @f$ \bar x_1 y_1 + \bar x_2 y_2 + \bar x_3 y_3 @f$ */
#define DOT3(x,y) (CONJ((x)[0]) * (y)[0] + CONJ((x)[1]) * (y)[1] + CONJ((x)[2]) * (y)[2])

/** @brief Compute squared 2-norm of a vector of  dimension 2, i.e.
 *   @f$ \lvert x \rvert^2 + \lvert y \rvert^2 @f$*/
#define NORMSQR2(x,y) (ABSSQR(x) + ABSSQR(y))

/** @brief Compute squared 2-norm of a vector of dimension 3, i.e.
 *   @f$ \lvert x \rvert^2 + \lvert y \rvert^2 + \lvert z \rvert^2 @f$*/
#define NORMSQR3(x,y,z) (ABSSQR(x) + ABSSQR(y) + ABSSQR(z))

/** @brief Compute 2-norm of a vector of dimension 2, i.e.
 *   @f$ \sqrt{\lvert x \rvert^2 + \lvert y \rvert^2} @f$*/
#define NORM2(x,y) REAL_SQRT(ABSSQR(x) + ABSSQR(y))

/** @brief Compute 2-norm of a vector of dimension 3, i.e.
 *   @f$ \sqrt{\lvert x \rvert^2 + \lvert y \rvert^2 + \lvert z \rvert^2} @f$*/
#define NORM3(x,y,z) REAL_SQRT(ABSSQR(x) + ABSSQR(y) + ABSSQR(z))

/** @brief Compute dot product of vectors of reals with dimension 2, i.e.
 *   @f$ x_1 y_1 + x_2 y_2 @f$ */
#define REAL_DOT2(x,y) ((x)[0] * (y)[0] + (x)[1] * (y)[1])

/** @brief Compute dot product of vectors of reals with dimension 3, i.e.
 *   @f$ x_1 y_1 + x_2 y_2 + x_3 y_3 @f$ */
#define REAL_DOT3(x,y) ((x)[0] * (y)[0] + (x)[1] * (y)[1] + (x)[2] * (y)[2])

/** @brief Compute squared 2-norm of a vector of reals with dimension 2, i.e.
 *   @f$ x^2 + y^2 @f$*/
#define REAL_NORMSQR2(x,y) (REAL_SQR(x) + REAL_SQR(y))

/** @brief Compute squared 2-norm of a vector of reals with dimension 3, i.e.
 *   @f$ x^2 + y^2 +z^2 @f$*/
#define REAL_NORMSQR3(x,y,z) (REAL_SQR(x) + REAL_SQR(y) + REAL_SQR(z))

/** @brief Compute 2-norm of a vector of reals with dimension 2, i.e.
 *   @f$ \sqrt{x^2 + y^2} @f$*/
#define REAL_NORM2(x,y) REAL_SQRT(REAL_SQR(x) + REAL_SQR(y))

/** @brief Compute 2-norm of a vector of reals with dimension 3, i.e.
 *   @f$ \sqrt{x^2 + y^2 + z^2} @f$*/
#define REAL_NORM3(x,y,z) REAL_SQRT(REAL_SQR(x) + REAL_SQR(y) + REAL_SQR(z))

/** @brief Compute the maximum @f$\max\{x,y\}@f$ of two real numbers @f$x@f$ and @f$y@f$. */
#define REAL_MAX(x, y) _h2_realmax(x, y)

/** @brief Compute the maximum @f$\max\{x,y,z\}@f$ of three real numbers @f$x,y@f$ and @f$z@f$. */
#define REAL_MAX3(x, y, z) _h2_realmax3(x, y, z)

/** @brief Compute the minimum @f$\min\{x,y\}@f$ of two real numbers @f$x@f$ and @f$y@f$. */
#define REAL_MIN(x, y) _h2_realmin(x, y)

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of three real numbers @f$x,y@f$ and @f$z@f$. */
#define REAL_MIN3(x, y, z) _h2_realmin3(x, y, z)

/* Macros for the type "uint" */

/** @brief Compute the maximum @f$\max\{x,y\}@f$ of two unsigned integers @f$x@f$ and @f$y@f$. */
#define UINT_MAX(x, y) _h2_uintmax(x, y)

/** @brief Compute the maximum @f$\max\{x,y,z\}@f$ of three unsigned integers @f$x,y@f$ and @f$z@f$. */
#define UINT_MAX3(x, y, z) _h2_uintmax3(x, y, z)

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of two unsigned integers @f$x,y@f$ and @f$z@f$. */
#define UINT_MIN(x, y) _h2_uintmin(x, y)

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of three unsigned integers @f$x,y@f$ and @f$z@f$. */
#define UINT_MIN3(x, y, z) _h2_uintmin3(x, y, z)

#ifdef __GNUC__
INLINE_PREFIX real _h2_rsqrt(real x) __attribute__((unused,const));
INLINE_PREFIX real _h2_abssqr(field x) __attribute__((unused,const));
INLINE_PREFIX field _h2_sgn1(field x) __attribute__((unused,const));
INLINE_PREFIX field _h2_fieldrand() __attribute__((unused,const));
INLINE_PREFIX real _h2_real_sqr(real x) __attribute__((unused,const));
INLINE_PREFIX real _h2_realmax(real a, real b) __attribute__((unused,const));
INLINE_PREFIX real _h2_realmax3(real a, real b, real c) __attribute__((unused,const));
INLINE_PREFIX real _h2_realmin(real a, real b) __attribute__((unused,const));
INLINE_PREFIX real _h2_realmin3(real a, real b, real c) __attribute__((unused,const));
INLINE_PREFIX real _h2_realrand() __attribute__((unused,const));
INLINE_PREFIX uint _h2_uintmax(uint a, uint b) __attribute__((unused,const));
INLINE_PREFIX uint _h2_uintmax3(uint a, uint b, uint c) __attribute__((unused,const));
INLINE_PREFIX uint _h2_uintmin(uint a, uint b) __attribute__((unused,const));
INLINE_PREFIX uint _h2_uintmin3(uint a, uint b, uint c) __attribute__((unused,const));
#else
INLINE_PREFIX real _h2_realrand();
#endif

/** @brief Compute the reciprocal square root @f$1.0/\sqrt{x}@f$ of a
 *  non-negative real nunber @f$x@f$.
 *  @param x Input value
 */
INLINE_PREFIX real _h2_rsqrt(real x) {
  return r_one / REAL_SQRT(x);
}

/** @brief Compute the square of the absolute value @f$|x|^2@f$ of a field element @f$x@f$.
 *
 *  @param x Real number @f$x@f$.
 *  @returns Square @f$|x|^2@f$. */
#ifdef USE_COMPLEX
INLINE_PREFIX real _h2_abssqr(field x) {
  real rx = REAL(x);
  real ix = IMAG(x);

  return rx * rx + ix * ix;
}
#else
INLINE_PREFIX real _h2_abssqr(field x) {
  return x * x;
}
#endif

/** @brief Compute the sign @f$\mathop{\rm sgn}(x)@f$ of a field
 *  element @f$x@f$ with the convention @f$\mathop{\rm sgn}(0)=1@f$.
 *
 *  @param x Field element @f$x@f$.
 *  @returns Sign of @f$x@f$, i.e., @f$x/|x|@f$ if @f$x\neq 0@f$ and @f$1@f$ if @f$x=0@f$. */
#ifdef USE_COMPLEX
INLINE_PREFIX field _h2_sgn1(field x) {
  real norm, rx, ix;

  if (ABSSQR(x) <= H2_ALMOST_ZERO) {
    return 1.0;
  } else {
    rx = REAL(x);
    ix = IMAG(x);
    norm = REAL_RSQRT(rx * rx + ix * ix);
    return rx * norm + ix * norm * I;
  }
}
#else
INLINE_PREFIX field _h2_sgn1(field x) {
  return (x < f_zero ? f_minusone : f_one);
}
#endif

/** @brief Compute a (pseudo-)random number.
 *
 *  @returns Random number with real and imaginary part in
 *    @f$[-1,1]@f$. */
INLINE_PREFIX field _h2_fieldrand() {
  field x;

  x = _h2_realrand();
#ifdef USE_COMPLEX
  x += f_i * _h2_realrand();
#endif

  return x;
}

/** @brief Compute the square @f$x^2@f$ of a real number @f$x@f$.
 *
 *  @param x Real number @f$x@f$.
 *  @returns Square @f$x^2@f$. */
INLINE_PREFIX real _h2_real_sqr(real x) {
  return x * x;
}

/** @brief Compute the maximum @f$\max\{x,y\}@f$ of two real numbers @f$x@f$ and @f$y@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @returns Maximum of @f$x@f$ and @f$y@f$. */
INLINE_PREFIX real _h2_realmax(real x, real y) {
  return (x < y ? y : x);
}

/** @brief Compute the maximum @f$\max\{x,y,z\}@f$ of three real numbers @f$x,y@f$ and @f$z@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @param z Third number.
 *  @returns Maximum of @f$x@f$, @f$y@f$, and @f$z@f$. */
INLINE_PREFIX real _h2_realmax3(real x, real y, real z) {
  return (x < y ? (y < z ? z : y) : (x < z ? z : x));
}

/** @brief Compute the minimum @f$\min\{x,y\}@f$ of two real numbers @f$x@f$ and @f$y@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @returns Minimum of @f$x@f$ and @f$y@f$. */
INLINE_PREFIX real _h2_realmin(real a, real b) {
  return (a < b ? a : b);
}

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of three real numbers @f$x,y@f$ and @f$z@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @param z Third number.
 *  @returns Minimum of @f$x@f$, @f$y@f$ and @f$z@f$. */
INLINE_PREFIX real _h2_realmin3(real x, real y, real z) {
  return (x < y ? (z < x ? z : x) : (z < y ? z : y));
}

/** @brief Compute a (pseudo-)random number.
 *
 *  @returns Random number in @f$[-1,1]@f$. */
INLINE_PREFIX real _h2_realrand() {
  return 2.0 * rand() / RAND_MAX - 1.0;
}

/** @brief Compute the maximum @f$\max\{x,y\}@f$ of two unsigned integers @f$x@f$ and @f$y@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @returns Maximum of @f$x@f$ and @f$y@f$. */
INLINE_PREFIX uint _h2_uintmax(uint x, uint y) {
  return (x < y ? y : x);
}

/** @brief Compute the maximum @f$\max\{x,y,z\}@f$ of three unsigned integers @f$x,y@f$ and @f$z@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @param z Third number.
 *  @returns Maximum of @f$x@f$, @f$y@f$ and @f$z@f$. */
INLINE_PREFIX uint _h2_uintmax3(uint x, uint y, uint z) {
  return (x < y ? (y < z ? z : y) : (x < z ? z : x));
}

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of two unsigned integers @f$x,y@f$ and @f$z@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @returns Minimum of @f$x@f$ and @f$y@f$. */
INLINE_PREFIX uint _h2_uintmin(uint x, uint y) {
  return (x < y ? x : y);
}

/** @brief Compute the minimum @f$\min\{x,y,z\}@f$ of three unsigned integers @f$x,y@f$ and @f$z@f$.
 *
 *  @param x First number.
 *  @param y Second number.
 *  @param z Third number.
 *  @returns Minimum of @f$x@f$, @f$y@f$ and @f$z@f$. */
INLINE_PREFIX uint _h2_uintmin3(uint x, uint y, uint z) {
  return (x < y ? (z < x ? z : x) : (z < y ? z : y));
}

/* ------------------------------------------------------------
 * Memory management
 * ------------------------------------------------------------ */

/** @brief Allocate heap storage.
 *
 *  @param sz Number of bytes.
 *  @returns Pointer to <tt>sz</tt> bytes. */
#define allocmem(sz) _h2_allocmem(sz,__FILE__,__LINE__)
/** @brief Allocate heap storage.
 *
 *  @param sz Number of bytes.
 *  @param filename Name of source file (used for error messages).
 *  @param line Line number in source file.
 *  @returns Pointer to <tt>sz</tt> bytes. */
HEADER_PREFIX void *
_h2_allocmem(size_t sz, const char *filename, int line);

/** @brief Allocate heap storage of type @ref uint.
 *
 *  @param sz Number of @ref uint variables.
 *  @returns Pointer to <tt>sz</tt> variables of type @ref uint. */
#define allocuint(sz) _h2_allocuint(sz,__FILE__,__LINE__)
/** @brief Allocate heap storage of type @ref uint.
 *
 *  @param sz Number of @ref uint variables.
 *  @param filename Name of source file (used for error messages).
 *  @param line Line number in source file.
 *  @returns Pointer to <tt>sz</tt> variables of type @ref uint. */
HEADER_PREFIX uint *
_h2_allocuint(size_t sz, const char *filename, int line);

/** @brief Allocate heap storage of type @ref real.
 *
 *  @param sz Number of @ref real variables
 *  @returns Pointer to <tt>sz</tt> variables of type @ref real. */
#define allocreal(sz) _h2_allocreal(sz,__FILE__,__LINE__)
/** @brief Allocate heap storage of type @ref real.
 *
 *  @param sz Number of @ref real variables.
 *  @param filename Name of source file (used for error messages).
 *  @param line Line number in source file.
 *  @returns Pointer to <tt>sz</tt> variables of type @ref real. */
HEADER_PREFIX real *
_h2_allocreal(size_t sz, const char *filename, int line);

/** @brief Allocate heap storage of type @ref field.
 *
 *  @param sz Number of @ref field variables.
 *  @returns Pointer to <tt>sz</tt> variables of type @ref field. */
#define allocfield(sz) _h2_allocfield(sz,__FILE__,__LINE__)
/** @brief Allocate heap storage of type <tt>field</tt>.
 *
 *  @param sz Number of @ref field variables.
 *  @param filename Name of source file (used for error messages).
 *  @param line Line number in source file.
 *  @returns Pointer to <tt>sz</tt> variables of type @ref field. */
HEADER_PREFIX field *
_h2_allocfield(size_t sz, const char *filename, int line);

/** @brief Allocate heap storage for a matrix.
 *
 *  \e Remark: For algebraic operations, consider @ref new_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Pointer to <tt>rows*cols</tt> variables of type @ref field. */
#define allocmatrix(rows,cols) _h2_allocmatrix(rows,cols,__FILE__,__LINE__)
/** @brief Allocate heap storage for a matrix with coefficients
 *  of type <tt>field</tt>.
 *
 *  \e Remark: For algebraic operations, consider using @ref new_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @param filename Name of source file (used for error messages).
 *  @param line Line number in source file.
 *  @returns Pointer to <tt>rows*cols</tt> variables of type @ref field. */
HEADER_PREFIX field *
_h2_allocmatrix(size_t rows, size_t cols, const char *filename, int line);

/** @brief Release allocated storage.
 *
 *  @param ptr Pointer to allocated storage. */
HEADER_PREFIX void
freemem(void *ptr);

/* ------------------------------------------------------------
 * Sorting
 * ------------------------------------------------------------ */

/** @brief Heapsort algorithm.
 *
 *  This function permutes the given array to put it into increasing
 *  order.
 *
 *  @param n Number of array elements.
 *  @param leq Return <tt>true</tt> if the array element corresponding
 *  to the first index is less than or equal to the one corresponding to
 *  the second index.
 *  @param swap Swap the array elements corresponding to the first and
 *  second index.
 *  @param data Array to be sorted. */
#define heapsort(n,leq,swap,data) _h2_heapsort(n,leq,swap,data)

/** @brief Heapsort algorithm.
 *
 *  This function permutes the given array to put it into increasing
 *  order.
 *
 *  @param n Number of array elements.
 *  @param leq Return <tt>true</tt> if the array element corresponding
 *  to the first index is less than or equal to the one corresponding to
 *  the second index.
 *  @param swap Swap the array elements corresponding to the first and
 *  second index.
 *  @param data Array to be sorted. */
HEADER_PREFIX void
_h2_heapsort(uint n, bool leq(uint, uint, void *),
    void swap(uint, uint, void *), void *data);

/* ------------------------------------------------------------
 * Timing
 * ------------------------------------------------------------ */

/** @brief Create a @ref stopwatch object.
 *
 *  @returns New @ref stopwatch object. */
HEADER_PREFIX pstopwatch
new_stopwatch();

/** @brief Delete a @ref stopwatch object.
 *
 *  @param sw Stopwatch object to be deleted. */
HEADER_PREFIX void
del_stopwatch(pstopwatch sw);

/** @brief Start a stopwatch.
 *
 *  This function stores the current time in a private
 *  variable that can be used subsequently to measure elapsed time.
 *
 *  @param sw Stopwatch to be started. */
HEADER_PREFIX void
start_stopwatch(pstopwatch sw);

/** @brief Stop a stopwatch.
 *
 *  This function stores the current time in a private
 *  variable that can be used subsequently to determine the time
 *  that has passed between calls to @ref start_stopwatch and
 *  @ref stop_stopwatch.
 *
 *  @param sw Stopwatch (@ref start_stopwatch has to have been called for this object at least once).
 *  @returns Elapsed time in seconds since stopwatch was started. */
HEADER_PREFIX real
stop_stopwatch(pstopwatch sw);

/* ------------------------------------------------------------
 *  Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
/** @brief Create a PDF canvas for Cairo drawing.
 *
 *  @remark The resulting <tt>cairo_t</tt> object should be destroyed
 *  using <tt>cairo_destroy</tt> in order to ensure that pending
 *  drawing operations are completed and the PDF file is closed.
 *
 *  @param filename Name of the PDF file.
 *  @param width Width of the canvas in pixels.
 *  @param height Height of the canvas in pixels.
 *  @returns <tt>cairo_t</tt> object that can be used in drawing operations. */
HEADER_PREFIX cairo_t *
new_cairopdf(const char *filename, double width, double height);

/** @brief Create a PNG canvas for Cairo drawing.
 *
 *  The resulting <tt>cairo_t</tt> object can be used in
 *  @ref write_cairopng to write the image to a PNG file.
 *
 *  @remark The <tt>cairo_t</tt> object should be destroyed
 *  using <tt>cairo_destroy</tt> once it is no longer needed.
 *
 *  @param width Width of the canvas in pixels.
 *  @param height Height of the canvas in pixels.
 *  @returns <tt>cairo_t</tt> object that can be used in drawing operations. */
HEADER_PREFIX cairo_t *
new_cairopng(uint width, uint height);

/** @brief Write Cairo image to PNG file.
 *
 *  @param cr Cairo context created by @ref new_cairopng.
 *  @param filename Name of the PNG file.
 *  @returns True if successful. */
HEADER_PREFIX bool
write_cairopng(cairo_t *cr, const char *filename);
#endif

/** @} */

/** @mainpage
 *
 *  The <tt>H2Lib</tt> package contains algorithms and data structures
 *  for working with hierarchical matrices [@cite HA99],
 *  [@cite GRHA02], [@cite HA09] and @f$\mathcal{H}^2@f$-matrices
 *  [@cite HAKHSA00], [@cite BOHA02], [@cite BO10].
 *  It is being developed in the Scientific Computing Group of
 *  Kiel University, since we require a software library that can
 *  be used both for teaching and research purposes and the existing
 *  libraries currently do not meet both requirements.
 *
 *  In order to offer a good basis for teaching courses on
 *  hierarchical matrices, the modules have been organized in a layered
 *  design that allows students to work with the lower layers (e.g.,
 *  for handling matrices and vectors) without having to worry about
 *  higher layers (e.g., approximative algebraic routines or sophisticated
 *  compression algorithms).
 *
 *  A researcher using <tt>H2Lib</tt> finds a relatively complete
 *  set of functions for creating and manipulating @f$\mathcal{H}@f$-
 *  and @f$\mathcal{H}^2@f$-matrices, e.g., for performing
 *  matrix-vector multiplications, approximative algebraic operations
 *  like multiplication, inversion and factorization, and functions
 *  for compressing and converting matrices between different
 *  representations.
 *
 *  For the sake of convenience, we have also included modules for
 *  a number of typical model problems, e.g., for boundary integral
 *  equations or elliptic partial differential equations, with the
 *  corresponding auxiliary modules for singular quadrature and
 *  simple grid management.
 *
 */

#endif
