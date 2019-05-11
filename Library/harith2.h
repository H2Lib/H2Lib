
/* ------------------------------------------------------------
   This is the file "harith2.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2015
   ------------------------------------------------------------ */

/** @file harith2.h
    @author Steffen B&ouml;rm
*/

#ifndef HARITH2_H
#define HARITH2_H

/** @addtogroup harith
 *  @{ */

/** @brief Accumulator for H-matrix products. */
typedef struct _haccum haccum;

/** @brief Pointer to @ref haccum object. */
typedef haccum *phaccum;

/** @brief Pointer to constant @ref haccum object. */
typedef const haccum *pchaccum;

/** @brief Unhandled H-matrix products. */
typedef struct _hprodentry hprodentry;

/** @brief Points to @ref hprodentry object. */
typedef hprodentry *phprodentry;

#include "harith.h"

/** @brief Accumulator for H-matrix products.
 *
 *  In order to reduce the number of SVD computations,
 *  we can accumulate updates to the same submatrix in an
 *  auxiliary structure.
 *  Low-rank updates are accumulated in an @ref rkmatrix object,
 *  while structured updates are entered into a list of
 *  @ref hprodentry objects. */
struct _haccum {
  /** @brief Target matrix. */
  phmatrix z;

  /** @brief Truncation mode. */
  pctruncmode tm;

  /** @brief Truncation accuracy. */
  real eps;
  
  /** @brief Accumulator for low-rank updates. */
  prkmatrix r;

  /** @brief Accumulator for structured updates. */
  phprodentry xy;
};

/** @brief Create a new H-matrix accumulator.
 *
 *  @param z Target matrix.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns Accumulator for <tt>z</tt> */
HEADER_PREFIX phaccum
new_haccum(phmatrix z, pctruncmode tm, real eps);

/** @brief Delete an accumulator.
 *
 *  @param ha Accumulator to be deleted. */
HEADER_PREFIX void
del_haccum(phaccum ha);

/** @brief Add a product @f$\alpha X Y@f$ to an accumulator.
 *
 *  @param alpha Scaling factor.
 *  @param xtrans Set if @f$X^*@f$ is to be used instead of @f$X@f$.
 *  @param x H-matrix @f$X@f$.
 *  @param ytrans Set if @f$Y^*@f$ is to be used instead of @f$Y@f$.
 *  @param y H-matrix @f$Y@f$.
 *  @param za Accumulator for the result. */
HEADER_PREFIX void
addproduct_haccum(field alpha,
		  bool xtrans, pchmatrix x, bool ytrans, pchmatrix y,
		  phaccum za);

/** @brief Create H-matrix accumulators for submatrices.
 *
 *  @param z Target matrix, row and column cluster have to conincide
 *     with <tt>ha->z->rc</tt> and <tt>ha->z->cc</tt>.
 *  @param za H-matrix accumulator.
 *  @returns Array containing accumulators for the submatrices of
 *     <tt>z</tt>, numbered as in <tt>z->son</tt>. */
HEADER_PREFIX phaccum *
split_haccum(phmatrix z, phaccum za);

/** @brief Flush an H-matrix accumulator, i.e., add all accumulated
 *    products to the target matrix and clear the accumulator.
 *
 *  @param ha H-matrix accumulator. */
HEADER_PREFIX void
flush_haccum(phaccum ha);

/** @brief Truncated addition of a matrix to a low-rank matrices
 *  in @ref rkmatrix representation,
 *  @f$B \gets \operatorname{trunc}(B+\alpha A,\epsilon)@f$.
 *
 *  @attention The matrix @f$A@f$ is overwritten by intermediate
 *  results.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix_destructive_rkmatrix(field alpha, bool atrans, pamatrix a,
         pctruncmode tm, real eps, prkmatrix b);

/** @brief Locally truncated addition of a matrix in @ref amatrix
 *  representation to a hierarchical matrix in @ref hmatrix representation,
 *  @f$B \gets \operatorname{blocktrunc}(B + \alpha A,\epsilon)@f$.
 *
 *  @attention The matrix @f$A@f$ is overwritten by intermediate
 *  results.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix_destructive_hmatrix(field alpha, bool atrans, pamatrix a,
        pctruncmode tm, real eps, phmatrix b);

/** @brief Multiply two H-matrices using accumulators,
 *  @f$Z \gets \operatorname{succtrunc}(Z + \alpha X Y,\epsilon)@f$.
 *
 *  Compared to the standard implementation @ref addmul_hmatrix, this
 *  algorithm uses H-matrix accumulators to reduce the number of
 *  singular value decompositions.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param xtrans Set if @f$X^*@f$ is to be used instead of @f$X@f$.
 *  @param x Hierarchical matrix @f$X@f$.
 *  @param ytrans Set if @f$Y^*@f$ is to be used instead of @f$Y@f$.
 *  @param y Hierarchical matrix @f$Y@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param z Target matrix @f$Z@f$. */
HEADER_PREFIX void
addmul2_hmatrix(field alpha,
		bool xtrans, pchmatrix x, bool ytrans, pchmatrix y,
		pctruncmode tm, real eps, phmatrix z);

/** @brief Solve a lower triangular system using accumulators,
 *  @f$X \gets \operatorname{succtrunc}(L^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(L^{-*} X, \epsilon)@f$.
 *
 *  @param aunit Set if @f$L@f$ has unit diagonal.
 *  @param atrans Set if @f$L^*@f$ instead of @f$L@f$ is to be used.
 *  @param a Matrix containing the lower triangular part of @f$L@f$.
 *     The remainder is not used.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xa Accumulator representation of @f$X@f$, will be overwritten
 *     by the solution. */
HEADER_PREFIX void
lowersolve_haccum(bool aunit, bool atrans, pchmatrix a,
		  bool xtrans, phaccum xa);

/** @brief Solve a lower triangular system using accumulators,
 *  @f$X \gets \operatorname{succtrunc}(L^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(L^{-*} X, \epsilon)@f$.
 *
 *  @param aunit Set if @f$L@f$ has unit diagonal.
 *  @param atrans Set if @f$L^*@f$ instead of @f$L@f$ is to be used.
 *  @param a Matrix containing the lower triangular part of @f$L@f$.
 *     The remainder is not used.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param x Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
lowersolve2_hmatrix_hmatrix(bool aunit, bool atrans, pchmatrix a,
			    pctruncmode tm, real eps,
			    bool xtrans, phmatrix x);

/** @brief Solve an upper triangular system using accumulators,
 *  @f$X \gets \operatorname{succtrunc}(R^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(R^{-*} X, \epsilon)@f$.
 *
 *  @param aunit Set if @f$R@f$ has unit diagonal.
 *  @param atrans Set if @f$R^*@f$ instead of @f$R@f$ is to be used.
 *  @param a Matrix containing the upper triangular part of @f$R@f$.
 *     The remainder is not used.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xa Accumulator representation of @f$X@f$, will be overwritten
 *     by the solution. */
HEADER_PREFIX void
uppersolve_haccum(bool aunit, bool atrans, pchmatrix a,
		  bool xtrans, phaccum xa);

/** @brief Solve an upper triangular system using accumulators,
 *  @f$X \gets \operatorname{succtrunc}(R^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(R^{-*} X, \epsilon)@f$.
 *
 *  @param aunit Set if @f$R@f$ has unit diagonal.
 *  @param atrans Set if @f$R^*@f$ instead of @f$R@f$ is to be used.
 *  @param a Matrix containing the lower triangular part of @f$R@f$.
 *     The remainder is not used.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param x Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
uppersolve2_hmatrix_hmatrix(bool aunit, bool atrans, pchmatrix a,
			    pctruncmode tm, real eps,
			    bool xtrans, phmatrix x);

/** @brief Compute the LR factorization using accumulators,
 *  @f$A \approx L R@f$.
 *
 *  The lower triangular part @f$L@f$ has unit diagonal, only its
 *  strict lower triangular part is stored in the strict lower
 *  triangular part of the source matrix.
 *
 *  The upper triangular part is stored in the upper triangular part
 *  of the source matrix.
 *
 *  @param aa Accumulator representation of @f$A@f$,
 *     will be overwritten with @f$L@f$ and @f$R@f$. */
HEADER_PREFIX void
lrdecomp_haccum(phaccum aa);

/** @brief Compute the LR factorization using accumulators,
 *  @f$A \approx L R@f$.
 *
 *  The lower triangular part @f$L@f$ has unit diagonal, only its
 *  strict lower triangular part is stored in the strict lower
 *  triangular part of the source matrix.
 *
 *  The upper triangular part is stored in the upper triangular part
 *  of the source matrix.
 *
 *  @param a Source matrix @f$A@f$,
 *     will be overwritten with @f$L@f$ and @f$R@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
lrdecomp2_hmatrix(phmatrix a, pctruncmode tm, real eps);

/** @brief Compute the Cholesky factorization using accumulators,
 *  @f$A \approx L L^*@f$.
 *
 *  The lower triangular part of @f$L@f$ is stored in the lower
 *  triangular part of the source matrix.
 *
 *  The strictly upper triangular part of the source matrix is
 *  not used.
 *
 *  @param aa Accumulator representation of @f$A@f$,
 *    lower triangular part will be overwritten */
HEADER_PREFIX void
choldecomp_haccum(phaccum aa);

/** @brief Compute the Cholesky factorization using accumulators,
 *  @f$A \approx L L^*@f$.
 *
 *  The lower triangular part of @f$L@f$ is stored in the lower
 *  triangular part of the source matrix.
 *
 *  The strictly upper triangular part of the source matrix is
 *  not used.
 *
 *  @param a Source matrix @f$A@f$, lower triangular part will be overwritten
 *    by @f$L@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
choldecomp2_hmatrix(phmatrix a, pctruncmode tm, real eps);

/** @} */

#endif
