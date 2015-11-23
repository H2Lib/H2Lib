/* ------------------------------------------------------------
 This is the file "truncation.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2014
 ------------------------------------------------------------ */

/** @file truncation.h
 *  @author Steffen B&ouml;rm */

#ifndef TRUNCATION_H
#define TRUNCATION_H

/** @defgroup truncation truncation
 *  @brief Auxiliary functions for truncation .
 *
 *  @{ */

/**
 * @brief Abbreviation for struct @ref _truncmode .
 */
typedef struct _truncmode truncmode;

/**
 * @brief Abbreviation for a pointer to a @ref truncmode object.
 */
typedef truncmode *ptruncmode;

/**
 * @brief Abbreviation for a pointer to a constant @ref truncmode object.
 */
typedef const truncmode *pctruncmode;

#include "realavector.h"
#include "settings.h"

/* ------------------------------------------------------------
 Truncation strategy
 ------------------------------------------------------------ */

/**
 * @brief Define different strategies used by various truncation and compression
 * algorithms for @ref hmatrix "hmatrices" and @ref h2matrix "h2matrices".
 */
struct _truncmode {
  /** @brief If set to <tt>true</tt> Frobenius instead of spectral norm is used. */
  bool frobenius;
  /** @brief If set to <tt>true</tt> Absolute instead of relative error is used. */
  bool absolute;
  /** @brief If set to <tt>true</tt> Blockwise error will be considered. */
  bool blocks;

  /** @brief Level-dependent tolerance factor. */
  real zeta_level;
  /** @brief Block-age-dependent tolerance factor */
  real zeta_age;
};

/* ------------------------------------------------------------
 Constructors for standard truncation strategies
 ------------------------------------------------------------ */

/**
 * @brief Create a new @ref truncmode object without any further initializations.
 *
 * @return A new, non specific @ref truncmode object is returned.
 */
HEADER_PREFIX ptruncmode
new_truncmode();

/**
 * @brief Free storage for @ref truncmode object.
 *
 * @param tm @ref truncmode "Truncmode" to be deleted.
 */
HEADER_PREFIX void
del_truncmode(ptruncmode tm);

/**
 * @brief Create a new @ref truncmode object with relative euclidean error norm.
 *
 * @return A new @ref truncmode object is returned using relative euclidean
 * error norm.
 */
HEADER_PREFIX ptruncmode
new_releucl_truncmode();

/**
 * @brief Create a new @ref truncmode object with relative frobenius error norm.
 *
 * @return A new @ref truncmode object is returned using relative frobenius
 * error norm.
 */
HEADER_PREFIX ptruncmode
new_relfrob_truncmode();

/**
 * @brief Create a new @ref truncmode object with block relative euclidean
 * error norm.
 *
 * @return A new @ref truncmode object is returned using block relative
 * euclidean error norm.
 */
HEADER_PREFIX ptruncmode
new_blockreleucl_truncmode();

/**
 * @brief Create a new @ref truncmode object with block relative frobenius
 * error norm.
 *
 * @return A new @ref truncmode object is returned using block relative
 * frobenius error norm.
 */
HEADER_PREFIX ptruncmode
new_blockrelfrob_truncmode();

/**
 * @brief Create a new @ref truncmode object with absolute euclidean error norm.
 *
 * @return A new @ref truncmode object is returned using absolute euclidean
 * error norm.
 */
HEADER_PREFIX ptruncmode
new_abseucl_truncmode();

/* ------------------------------------------------------------
 Find minimal acceptable rank
 ------------------------------------------------------------ */

/**
 * @brief Determine new rank for truncation routines.
 *
 * Depending on the truncation mode given by <tt>tm</tt>, a truncation accuracy
 * given by <tt>eps</tt> and the singular values given by the @ref avector
 * <tt>sigma</tt> the function returns the new suitable rank.
 *
 * @param tm Truncation strategy
 * @param eps Trucation accuracy
 * @param sigma @ref realavector "Realavector" containing the singular values.
 * @return Returns the new suitable rank.
 */
HEADER_PREFIX uint
findrank_truncmode(pctruncmode tm, real eps, pcrealavector sigma);

/**
 * @}
 */
#endif
