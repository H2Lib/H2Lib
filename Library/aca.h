/* ------------------------------------------------------------
 This is the file "aca.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2012
 ------------------------------------------------------------ */

/**
 * @file aca.h
 * @author Sven Christophersen
 * @date 2012
 */

#ifndef ACA_H_
#define ACA_H_

/* C STD LIBRARY */

/* CORE 0 */
#include "settings.h"
/* CORE 1 */
#include "amatrix.h"
/* CORE 2 */
#include "rkmatrix.h"
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/** \defgroup aca aca
 *  @brief This modules provides different kind of adaptive cross approximation
 *  techniques.
 *  @{ */

/** @brief Matrix entry callback.
 *
 *  Used to evaluate submatrices of an implicitly given matrix.
 *
 *  A typical application would be the computation of entries of
 *  a boundary element matrix by appropriate quadrature rules.
 *
 *  @param ridx Array of row indices. Should have at least <tt>N->rows</tt>
 *    elements if <tt>ntrans</tt> is not set and <tt>N->cols</tt> elements if
 *    it is set.
 *  @param cidx Array of column indices. Should have  at least <tt>N->cols</tt>
 *    elements if <tt>ntrans</tt> is not set and <tt>N->rows</tt> elemenrs if
 *    it is set.
 *  @param data Arbitrary data the callback function might require
 *    to complete its task, e.g., geometry data.
 *  @param ntrans Set if the transposed matrix is to be filled.
 *  @param N Target matrix. */
typedef void (*matrixentry_t)(const uint *ridx, const uint *cidx, void *data,
    const bool ntrans, pamatrix N);

/**
 * @brief This routine computes the adaptive cross approximation using full
 * pivoting of a given matrix @f$ A @f$.
 *
 * The matrix @f$A@f$ will be overwritten in part by the Schur complement,
 * while the matrix @f$R@f$ will be filled with the low-rank approximation
 * @f$R=C D^*@f$ of the original @f$R@f$.
 * After completion, we have
 * @f[
 * A \approx C \, D^*
 * @f]
 * with a guaranteed error of
 * @f[
 * \lVert A - C \, D^* \rVert \leq \epsilon \, \lVert A \rVert \, .
 * @f]
 *
 * @param A Input matrix for which the ACA has to be applied to.
 * @param accur Accuracy defining how good the approximation has to be relative
 * to the input matrix.
 * @param ridx Returns an array of row pivot indices, if <tt>ridx != NULL</tt>.
 * @param cidx Returns an array of column pivot indices, if <tt>cidx != NULL</tt>.
 * @param R The resulting low rank matrix is returned via <tt>R</tt> .
 */
HEADER_PREFIX void decomp_fullaca_rkmatrix(pamatrix A, const real accur,
    uint **ridx, uint **cidx, prkmatrix R);

/**
 * @brief This routine computes the adaptive cross approximation using partial
 * pivoting of an implicitly given matrix @f$ A @f$.
 *
 * After applying the ACA-algorithm we have
 * @f[
 * A \approx C \, D^*
 * @f]
 * with an estimated error of
 * @f[
 * \lVert A - C \, D^* \rVert \leq \epsilon \, \lVert A \rVert \, .
 * @f]
 *
 * @param entry This callback function implicitly defines the matrix @f$ A @f$.
 * It take two arrays of row and column indices, a void-pointer to some data-object
 * needed for computing matrix entries and a flag determined whether we want
 * entries of the original or of the transposed matrix @f$ A @f$. Length of <tt>
 * ridx</tt> and <tt>cidx</tt> is determined by the rows and columns of matrix
 * <tt>A</tt>.
 * @param data An additional void-pointer to some data-object that will be needed
 * by <tt>entry</tt> to compute the matrix entries.
 * @param ridx An array of all row indices defining the complete matrix @f$ A @f$
 * as in @ref decomp_fullaca_rkmatrix.
 * @param rows Number of rows for the implicit matrix and therefore the length of
 * <tt>ridx</tt>.
 * @param cidx An array of all column indices defining the complete matrix @f$ A @f$
 * as in @ref decomp_fullaca_rkmatrix.
 * @param cols Number of columns for the implicit matrix and therefore the length of
 * <tt>cidx</tt>.
 * @param accur Accuracy defining how good the approximation has to be relative
 * to the input matrix.
 * @param rpivot Returns an array of row pivot indices.
 * @param cpivot Returns an array of column pivot indices.
 * @param R The resulting low rank matrix is returned via <tt>R</tt> .
 */
HEADER_PREFIX void
decomp_partialaca_rkmatrix(matrixentry_t entry, void *data, const uint *ridx,
    const uint rows, const uint *cidx, const uint cols, real accur,
    uint **rpivot, uint **cpivot, prkmatrix R);

/**
 * @brief Copies the lower triangular part of a matrix <tt>A</tt> to a matrix <tt>B</tt>
 * after applying the row pivoting denoted by <tt>xi</tt>.
 *
 * for @f$ A \in \mathbb R^{\ell \times m}@f$ and
 * @f$ B \in \mathbb R^{n \times m}@f$ with @f$ \ell \geq n@f$ it holds
 *
 * @f[
 * B_{i,j} = A_{\xi(i),j} \quad \text{for all } i \in \{1,\ldots,n \}, \
 *    j \in \{ 1,\ldots,i \}.
 * @f]
 *
 * @param unit If set to true the diagonal will be set to 1.0 instead of the
 *        corresponding entry of <tt>A</tt>.
 * @param A Source matrix.
 * @param xi Array denoting the row permutation.
 * @param B Target matrix containing the lower triangular part of <tt>A</tt> with
 * respect to <tt>xi</tt>.
 *
 * @attention The rows of <tt>B</tt> has to match the length of the array
 * <tt>xi</tt>.
 */
HEADER_PREFIX void
copy_lower_aca_amatrix(bool unit, pcamatrix A, uint *xi, pamatrix B);

/**
 * @brief Copies the upper triangular part of a matrix <tt>A</tt> to a matrix <tt>B</tt>
 * after applying the row pivoting denoted by <tt>xi</tt>.
 *
 * for @f$ A \in \mathbb R^{\ell \times m}@f$ and
 * @f$ B \in \mathbb R^{n \times m}@f$ with @f$ \ell \geq n@f$ it holds
 *
 * @f[
 * B_{i,j} = A_{\xi(i),j} \quad \text{for all } i \in \{1,\ldots,n \}, \
 *    j \in \{ i,\ldots,m \}.
 * @f]
 *
 * @param unit If set to true the diagonal will be set to 1.0 instead of the
 *        corresponding entry of <tt>A</tt>.
 * @param A Source matrix.
 * @param xi Array denoting the row permutation.
 * @param B Target matrix containing the upper triangular part of <tt>A</tt> with
 * respect to <tt>xi</tt>.
 *
 * @attention The rows of <tt>B</tt> has to match the length of the array
 * <tt>xi</tt>.
 */
HEADER_PREFIX void
copy_upper_aca_amatrix(bool unit, pcamatrix A, uint *xi, pamatrix B);

/** @} */

#endif /* ACA_H_ */

