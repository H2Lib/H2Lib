
/* ------------------------------------------------------------
 * This is the file "settings.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file settings.h
 *  @author Steffen B&ouml;rm
 */

#ifndef SETTINGS_H
#define SETTINGS_H

/** @defgroup settings settings
 *  @brief Fundamental types and macros.
 *  @{ */

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#ifdef USE_COMPLEX
#include <complex.h>
#endif

#ifdef USE_SIMD
#include "simd.h"
#endif

/* ------------------------------------------------------------
 * Compilation settings
 * ------------------------------------------------------------ */

/** @brief Prefix for inline functions. */
#ifdef __cplusplus
#define INLINE_PREFIX inline
#else
#define INLINE_PREFIX static
#endif

/** @brief Prefix for function declarations. */
#ifdef __cplusplus
#define HEADER_PREFIX extern "C"
#else
#define HEADER_PREFIX
#endif

/** @brief Prefix for external declarations. */
#ifdef __cplusplus
#define IMPORT_PREFIX extern "C"
#else
#define IMPORT_PREFIX
#endif

/* ------------------------------------------------------------
 * Types
 * ------------------------------------------------------------ */

/** @brief Unsigned integer type.
 *
 *  This type is mostly used to access components of arrays,
 *  vectors and matrices. */
typedef unsigned uint;

/** @brief Signed integer constant zero. */
extern const int i_zero;

/** @brief Signed integer constant one. */
extern const int i_one;

/** @brief Unsigned integer constant zero. */
extern const uint u_zero;

/** @brief Unsigned integer constant one. */
extern const uint u_one;

/** @brief Unsigned long type.
 *
 *  This type is used to access components of particularly large
 *  arrays, e.g., matrices in column-major array representation. */
typedef size_t longindex;

/** @brief @ref real floating point type.
 *
 *  This type is used, e.g., for geometric coordinates, norms
 *  and diagonal elements of self-adjoint matrices. */
#ifdef USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

/**
 * @brief Fallback vector length for float if vectorization is not enabled.
 */
#ifndef VFLOAT
#define VFLOAT 1
#endif

/**
 * @brief Fallback vector length for double if vectorization is not enabled.
 */
#ifndef VDOUBLE
#define VDOUBLE 1
#endif

/**
 * @brief Fallback vector length for real if vectorization is not enabled.
 */
#ifndef VREAL
#define VREAL 1
#endif

#ifdef USE_FLOAT
/** Relative tolerance for run-time checks. */
#define H2_CHECK_TOLERANCE 1.0e-6

/** Bound for determining when a number is essentially zero. */
#define H2_ALMOST_ZERO 1e-30
#else
/** Relative tolerance for run-time checks. */
#define H2_CHECK_TOLERANCE 1.0e-12

/** Bound for determining when a number is essentially zero. */
#define H2_ALMOST_ZERO 1e-300
#endif

/**
 * @brief Prefix that is needed when reading numbers via 'scanf' function.
 */
#ifdef USE_FLOAT
#define SCANF_PREFIX ""
#else
#define SCANF_PREFIX "l"
#endif

/** @brief Pointer to @ref real array. */
typedef real *preal;

/** @brief Pointer to constant @ref real array. */
typedef const real *pcreal;

/** @brief @ref real constant zero */
extern const real r_zero;

/** @brief @ref real constant one. */
extern const real r_one;

/** @brief @ref real constant minus one. */
extern const real r_minusone;

/** @brief Field type.
 *
 *  This type is used in the linear algebra modules to represent
 *  the coefficients of matrices and vectors. */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
typedef float _Complex field;
#else
typedef float field;
#endif
#else
#ifdef USE_COMPLEX
typedef double _Complex field;
#else
typedef double field;
#endif
#endif

/**
 * @brief Fallback vector length for field if vectorization is not enabled.
 */
#ifndef VFIELD
#define VFIELD 1
#endif

/** @brief Pointer to @ref field array. */
typedef field *pfield;

/** @brief Pointer to constant @ref field array. */
typedef const field *pcfield;

/** @brief @ref field constant zero */
extern const field f_zero;

/** @brief @ref field constant one. */
extern const field f_one;

/** @brief @ref field constant minus one. */
extern const field f_minusone;

#ifdef USE_COMPLEX
/** @brief @ref field constant for the imaginary number. */
extern const field f_i;
#endif

#ifndef USE_COMPLEX
#ifdef I
#undef I
#endif
#define I 0.0
#endif

/** @brief String constant that expresses that a matrix should @b not be
 *   transposed. */
extern const char *_h2_ntrans;

/** @brief String constant that expresses that a matrix should be transposed. */
extern const char *_h2_trans;

/** @brief String constant that expresses that a matrix should be conjugated
 * and transposed. In case @ref field equals @ref real this is equivalent to
 * @ref _h2_trans. */
extern const char *_h2_adj;

/** @brief String constant that expresses that a matrix should be applied from
 *  the left side of another. */
extern const char *_h2_left;

/** @brief String constant that expresses that a matrix should be applied from
 *  the right side of another. */
extern const char *_h2_right;

/** @brief String constant that expresses that only the lower part of a matrix
 * should be used. */
extern const char *_h2_lower;

/** @brief String constant that expresses that only the upper part of a matrix
 * should be used. */
extern const char *_h2_upper;

/** @brief String constant that expresses that a matrix has implicitly given
 * unit diagonal entries. */
extern const char *_h2_unit;

/** @brief String constant that expresses that a matrix doesn't have implicitly
 * given  unit diagonal entries. */
extern const char *_h2_nonunit;

/** @brief String constant that expresses that all eigenvectors should be
 * computed by some eigenvalue/eigenvectors routine. */
extern const char *_h2_vectors;

/** @brief String constant that expresses that only the skinny eigenvectors
 * should be computed by some eigenvalue/eigenvectors routine. */
extern const char *_h2_skinnyvectors;

/** @brief String constant that expresses that no eigenvectors should be
 * computed by some eigenvalue/eigenvectors routine. */
extern const char *_h2_novectors;

/****************************************************
 * Define conversion specifier and argument macros
 * for different @ref field types.
 ****************************************************/

/**
 * @brief Macro that simplifies the definition of correct conversion specifier
 * for @ref field type values used in 'printf' function.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define FIELD_CS(format, cs) "%" #format #cs " + %" #format #cs "i"
#else
#define FIELD_CS(format, cs) "%" #format #cs
#endif
#else
#ifdef USE_COMPLEX
#define FIELD_CS(format, cs) "%" #format #cs " + %" #format #cs "i"
#else
#define FIELD_CS(format, cs) "%" #format #cs
#endif
#endif

/**
 * @brief Macro that simplifies the definition of correct conversion specifier
 * for @ref field type values used in 'scanf' function.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define FIELD_SCANF_CS(format, cs) "%" #format #cs " + %" #format #cs "i"
#else
#define FIELD_SCANF_CS(format, cs) "%" #format #cs
#endif
#else
#ifdef USE_COMPLEX
#define FIELD_SCANF_CS(format, cs) "%l" #format #cs " + %l" #format #cs "i"
#else
#define FIELD_SCANF_CS(format, cs) "%l" #format #cs
#endif
#endif

/**
 * @brief Macro that simplifies argument passing for 'printf' like functions.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define FIELD_ARG(z) crealf(z), cimagf(z)
#else
#define FIELD_ARG(z) z
#endif
#else
#ifdef USE_COMPLEX
#define FIELD_ARG(z) creal(z), cimag(z)
#else
#define FIELD_ARG(z) z
#endif
#endif

/**
 * @brief Macro that simplifies argument passing for 'printf' like functions
 * when the address of the parameters is needed.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define FIELD_ADDR(z) (real*)(z), ((real*)(z)+1)
#else
#define FIELD_ADDR(z) (real*)(z)
#endif
#else
#ifdef USE_COMPLEX
#define FIELD_ADDR(z) (real*)(z), ((real*)(z)+1)
#else
#define FIELD_ADDR(z) (real*)(z)
#endif
#endif

/** @brief All possible types of matrices.
 *
 *  Passing a flag of this type makes function using generic matrix types more
 *  legible.
 */
typedef enum {
  /** @brief Enum value representing an @ref _amatrix "amatrix". */
  AMATRIX = 0, //!< AMATRIX
  /** @brief Enum value representing an @ref _hmatrix "hmatrix". */
  HMATRIX = 1, //!< HMATRIX
  /** @brief Enum value representing an @ref _h2matrix "h2matrix". */
  H2MATRIX = 2, //!< H2MATRIX
  /** @brief Enum value representing an @ref _sparsematrix "sparsematrix". */
  SPARSEMATRIX = 3, //!< SPARSEMATRIX
  /** @brief Enum value representing an @ref _dh2matrix "dh2matrix". */
  DH2MATRIX = 4 //!< DH2MATRIX
} matrixtype;

/** @} */

#endif
