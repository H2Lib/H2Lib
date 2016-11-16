
/* ------------------------------------------------------------
 * This is the file "dh2compression.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file dh2compression.h
 *  @author Steffen B&ouml;rm
 */

#ifndef DH2COMPRESSION_H
#define DH2COMPRESSION_H

/** @defgroup dh2compression dh2compression
 *  @brief Functions for turning an @ref amatrix into an @ref dh2matrix.
 *
 *  The functions in this module can be used to convert array
 *  matrices into directional @f$\mathcal{H}^2@f$-matrices.
 *  @{ */

#include "dh2matrix.h"
#include "truncation.h"

/* ------------------------------------------------------------
 * Compute adaptive cluster bases for dense matrices
 * ------------------------------------------------------------ */

/** @brief Construct a directional row basis for an array matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param b Directional block tree.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New row cluster basis. */
pdclusterbasis
buildrow_amatrix_dclusterbasis(pcamatrix G, pcdblock b,
			       pctruncmode tm, real eps);

/** @brief Construct a directional column basis for an array matrix.
 *
 *  @param G Original matrix @f$G@f$.
 *  @param b Directional block tree.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns New column cluster basis. */
pdclusterbasis
buildcol_amatrix_dclusterbasis(pcamatrix G, pcdblock b,
			       pctruncmode tm, real eps);

/* ------------------------------------------------------------
 * Approximate dense matrix in new cluster bases
 * ------------------------------------------------------------ */

/** @brief Compute the optimal coupling matrix for a given matrix
 *    and given orthogonal row and column cluster bases.
 *
 *  @param a Source matrix, not in cluster ordering.
 *  @param rb Orthogonal row cluster basis.
 *  @param rd Direction for row basis.
 *  @param cb Orthogonal column cluster basis.
 *  @param cd Direction for column basis.
 *  @param s Will be filled with optimal coefficients to approximate
 *    the given matrix. */
void
collectdense_dh2matrix(pcamatrix a,
		       pcdclusterbasis rb, uint rd,
		       pcdclusterbasis cb, uint cd,
		       pamatrix s);

/** @brief Fill a given @ref dh2matrix object with the optimal
 *  coefficients for the approximation of a given matrix.
 *
 *  @param a Source matrix, not in cluster ordering.
 *  @param h2 Target matrix with orthogonal row and column cluster
 *    basis, nearfield and coupling matrices will be filled with
 *    coefficients approximation <tt>a</tt>. */
void
projectdense_dh2matrix(pcamatrix a, pdh2matrix h2);

/** @} */

#endif

