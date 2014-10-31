/* ------------------------------------------------------------
 This is the file "settings.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

/** @file settings.h
 @author Steffen B&ouml;rm
 */

#ifndef SETTINGS_H
#define SETTINGS_H

/** @defgroup settings settings
 *  @brief Fundamental types and macros.
 *  @{ */

#include <math.h>

/* ------------------------------------------------------------
 Compilation settings
 ------------------------------------------------------------ */

/** @brief Prefix for inline functions. */
#ifdef __cplusplus
#define INLINE_PREFIX static extern "C"
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
 Types
 ------------------------------------------------------------ */

/** @brief Boolean type. */
typedef unsigned short bool;

/** @brief Boolean constant <tt>true</tt>. */
extern const bool true;

/** @brief Boolean constant <tt>false</tt>. */
extern const bool false;

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
typedef unsigned longindex;

/** @brief @ref real floating point type.
 *
 *  This type is used, e.g., for geometric coordinates, norms
 *  and diagonal elements of self-adjoint matrices. */
typedef double real;

/** @brief Pointer to @ref real array. */
typedef real *preal;

/** @brief Pointer to constant @ref real array. */
typedef const real *pcreal;

/** @brief Field type.
 *
 *  This type is used in the linear algebra modules to represent
 *  the coefficients of matrices and vectors. */
typedef double field;

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

/** @brief All possible types of matrices.
 *
 *  Passing a flag of this type makes function using generic matrix type more
 *  legible.
 */

typedef enum {
/** @brief Enum value representing an @ref _amatrix "amatrix". */
 AMATRIX = 0,
 /** @brief Enum value representing an @ref _hmatrix "hmatrix". */
 HMATRIX = 1,
 /** @brief Enum value representing an @ref _h2matrix "h2matrix". */
 H2MATRIX = 2
 } matrixtype;

 /** @} */

#endif
