
/* ------------------------------------------------------------
   This is the file "h2update.h" of the H2Lib package.
   All rights reserved, Knut Reimer 2010
   ------------------------------------------------------------ */

/** @file h2update.h
 *  @author Knut Reimer
 */

#ifndef H2UPDATE_H
#define H2UPDATE_H

/** @defgroup h2update h2update
 *  @brief Functions for computing a low rank update for an @ref h2matrix.
 *
 *  The functions in this module can be used to compute low rank updates 
 *  for @f$\mathcal{H}^2@f$-matrices.
 *  The module also provides functions to initialise the auxilliary structures
 *  required by the low rank update and the arithmetic functions in @ref h2arith.
 *  @{ */

#include "rkmatrix.h"
#include "clusteroperator.h"
#include "h2matrix.h"
#include "h2compression.h"


/**
 *  @brief Computes the Euclidean norm of the extended coupling matrix of the low
 *  rank update multiplied by weight matrices.
 *
 *  Depending on given row and column @ref clusterbasis of the @ref uniform matrix
 *  this function computes one of the following four values via
 *  vector iteration:
 *  @f[
 *  \left\lVert \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} \right\rVert_2, \quad
 *  \left\lVert C_t \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} \right\rVert_2, \quad
 *  \left\lVert \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} C_s^* \right\rVert_2 \quad
 *  \text{or} \quad \left\lVert C_t \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} C_s^* \right\rVert_2
 *  @f]
 *  with @f$C_t@f$<tt>=u->rb->Z</tt> and @f$C_s@f$<tt>=u->cb->Z</tt>.
 *  @f$C_t@f$ is omitted if <tt>u->rb->Z==NULL</tt> and @f$C_s@f$ is omitted if <tt>u->cb->Z==NULL</tt>
 *
 *  @param u @ref _uniform "uniform" object
 *  @param k rank of the low rank matrix
 * 
 *  @return Approximation of the Euclidean norm of one of the matrices
 *    mentioned above.
 */
HEADER_PREFIX real
norm2_rkupdate_uniform(puniform u, uint k);

/**
 *  @brief Computes the Frobenius norm of the extended coupling matrix of the low
 *  rank update multiplied by weight matrices.
 *
 *   Depending on given row and column @ref clusterbasis of the @ref uniform matrix
 *  this functions computed one of the following four values:
 *   @f[
 *  \left\lVert \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} \right\rVert_F, \quad
 *  \left\lVert C_t \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} \right\rVert_F, \quad
 *  \left\lVert \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} C_s^* \right\rVert_F \quad
 *  \text{or} \quad \left\lVert C_t \begin{pmatrix} S_b & \\ & I_k \end{pmatrix} C_s^* \right\rVert_F
 *   @f]
 *  with @f$C_t@f$<tt>=u->rb->Z</tt> and @f$C_s@f$<tt>=u->cb->Z</tt>.
 *  @f$C_t@f$ is omitted if <tt>u->rb->Z==NULL</tt> and @f$C_s@f$ is omitted if <tt>u->cb->Z==NULL</tt>
 *
 *  @param u @ref _uniform "uniform" object
 *  @param k rank of the low rank matrix
 * 
 *  @return Approximation of the Frobenius norm of one of the matrices
 *    mentioned above.
 */
HEADER_PREFIX real
normfrob_rkupdate_uniform(puniform u, uint k);

/**
 *  @brief Computes the low rank update @f$ G \gets G + R @f$
 * 
 *  @param R Low-rank matrix @f$R@f$.
 *  @param Gh2 Target matrix @f$G@f$.
 *  @param rwf has to be the father of the total weights of the row clusterbasis of C,\n
 *    e.g. initialised by prepare_row_clusteroperator
 *  @param cwf has to be the father of the total weights of the col clusterbasis of C,\n
 *    e.g. initialised by prepare_col_clusteroperator
 *  @param tm options of truncation
 *  @param eps tolerance of truncation
 */
HEADER_PREFIX void
rkupdate_h2matrix(prkmatrix R, ph2matrix Gh2, pclusteroperator rwf, pclusteroperator cwf,
                   ptruncmode tm, real eps);

/**
 * @brief Prepares the weights of the row clusterbasis used by @ref rkupdate_h2matrix and the arithmetic functions in @ref h2arith  
 * 
 *  @param rb row clusterbasis
 *  @param cb column clusterbasis
 *  @param tm options for the total weights
 *  @return new clusteroperator
 */
HEADER_PREFIX pclusteroperator
prepare_row_clusteroperator(pclusterbasis rb, pclusterbasis cb, ptruncmode tm);

/** @brief Prepares the weights of the column clusterbasis used by @ref rkupdate_h2matrix and the arithmetic functions in @ref h2arith  
 * 
 *  @param rb row clusterbasis
 *  @param cb column clusterbasis
 *  @param tm options for the total weights
 *  @return new clusteroperator
 */
HEADER_PREFIX pclusteroperator
prepare_col_clusteroperator(pclusterbasis rb, pclusterbasis cb, ptruncmode tm);

/** @} */

#endif
