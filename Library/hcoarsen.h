/* ------------------------------------------------------------
 This is the file "hcoarsen.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2014
 ------------------------------------------------------------ */

/**
 * @file hcoarsen.h
 * @author Sven Christophersen
 */

#ifndef HCOARSEN_H
#define HCOARSEN_H

#include "harith.h"

/** @defgroup hcoarsen hcoarsen
 *  @brief Coarsening of hierarchical matrices.
 *  @{ */

/**
 * @brief Coarsen the block structure of a @ref hmatrix.
 *
 * For a @ref hmatrix consisting only of admissible or inadmissible leafs
 * it can be advantageous to store the @ref hmatrix as a single low rank matrix
 * with a slightly higher rank.
 *
 * If we consider a 2x2 block matrix
 * @f[
 * G = \begin{pmatrix} A_{1} B_{1}^* & A_{2} B_{2}^* \\
 *                     A_{3} B_{3}^* & A_{4} B_{4}^*
 *     \end{pmatrix}
 * @f]
 * with a low rank representation @f$ A_{i} B_{i}^*, i \in \{ 1,\ldots,4 \}@f$
 * with each having a rank of at most @f$ k @f$.
 * (Inadmissible blocks can be written as low rank matrices with
 * @f$ G_{i} = G_{i} I = A_{i} B_{i}^*@f$ or
 * @f$ G_{i} = I G_{i} = A_{i} B_{i}^*@f$.)
 *
 * We can rewrite this as a rank @f$ 4 k @f$ low rank matrix
 * @f[
 *  \widehat G = \begin{pmatrix} A_{1} & A_{2} \\
 *                     & & A_{3} & A_{4}
 *      \end{pmatrix}
 *      \begin{pmatrix} B_{1}^* &  \\
 *                      & B_{2}^*  \\
 *                      B_{3}^* &  \\
 *                      & B_{4}^*
 *      \end{pmatrix}
 * @f]
 *
 * Now an approximation of @f$ \widehat G@f$ via SVD using @ref trunc_rkmatrix
 * is computed:
 * @f[ \widetilde G = \operatorname{trunc}(\widehat G)
 *                  = \widetilde A \widetilde B^*@f].
 *
 * If the size of @f$ \widetilde G@f$ is smaller than the size of @f$ G @f$,
 * then @f$ G @f$ will be replaced by @f$ \widetilde G @f$. Otherwise
 * @f$ \widetilde G@f$ is rejected.
 *
 * If <tt>recursive == true</tt> holds, this process is repeated for father blocks
 * as long as they only consist of leaf blocks aswell.
 *
 * @param G Input @ref hmatrix. Will be changed during the coarsening process.
 * @param tm Truncation mode.
 * @param eps Accuracy for low rank truncation.
 * @param recursive Flag to indicate whether the coarsening algorithm should
 * be applied to the son blocks aswell or not.
 */
HEADER_PREFIX void
coarsen_hmatrix(phmatrix G, ptruncmode tm, real eps, bool recursive);

/**
 * @}
 */

#endif
