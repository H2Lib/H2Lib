/* ------------------------------------------------------------
 * This is the file "matrixnorms.h" of the H2Lib package.
 * All rights reserved, Sven Christophersen 2016
 * ------------------------------------------------------------ */

/** @file matrixnorms.h
 *  @author Sven Christophersen
 */

#ifndef MATRIXNORMS_H_
#define MATRIXNORMS_H_

#include "basic.h"
#include "krylov.h"
#include "amatrix.h"
#include "rkmatrix.h"
#include "hmatrix.h"
#include "h2matrix.h"
#include "dh2matrix.h"
#include "sparsematrix.h"
#include "factorizations.h"
#include "harith.h"

/** @defgroup matrixnorms matrixnorms
 *  @brief Estimators for spectral norm of difference form differnt types of
 *         matrices @f$A@f$ and @f$B@f$.
 *         Approximation of the form @f$\lVert A-B \rVert_2@f$ or
 *         @f$\lVert A-B \rVert_2@f$ with some factorized matrix @f$B@f$ or
 *         @f$\lVert I - B^{-1} A \rVert_2@f$ with some facotized matrix @f$B@f$
 *         can be computed.
 *
 *  @{ */

/****************************************************
 * Norm2diff for amatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Low rank matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix_rkmatrix(pcrkmatrix a, pcamatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Sparse matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix_sparsematrix(pcsparsematrix a, pcamatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Hierarchical matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix_hmatrix(pchmatrix a, pcamatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix_h2matrix(pch2matrix a, pcamatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix_dh2matrix(pcdh2matrix a, pcamatrix b);

/****************************************************
 * Norm2diff for rkmatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Sparse matrix @f$A@f$.
 *  @param b Low rank matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_rkmatrix_sparsematrix(pcsparsematrix a, pcrkmatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Hierarchical matrix @f$A@f$.
 *  @param b Low rank matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_rkmatrix_hmatrix(pchmatrix a, pcrkmatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Low rank matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_rkmatrix_h2matrix(pch2matrix a, pcrkmatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Low rank matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_rkmatrix_dh2matrix(pcdh2matrix a, pcrkmatrix b);

/****************************************************
 * Norm2diff for sparsematrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Hierarchical matrix @f$A@f$.
 *  @param b Sparse matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_sparsematrix_hmatrix(pchmatrix a, pcsparsematrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Sparse matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_sparsematrix_h2matrix(pch2matrix a, pcsparsematrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Sparse matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_sparsematrix_dh2matrix(pcdh2matrix a, pcsparsematrix b);

/****************************************************
 * Norm2diff for hmatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Hierarchical matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_hmatrix_h2matrix(pch2matrix a, pchmatrix b);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b Hierarchical matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_hmatrix_dh2matrix(pcdh2matrix a, pchmatrix b);

/****************************************************
 * Norm2diff for h2matrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b @f$\mathcal H^2@f$ matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_h2matrix_dh2matrix(pcdh2matrix a, pch2matrix b);

/****************************************************
 * Norm2diff_pre_matrix for amatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_amatrix(pcamatrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_amatrix_hmatrix(pcamatrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_amatrix(pcamatrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_amatrix_hmatrix(pcamatrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_pre_matrix for hmatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_hmatrix_amatrix(pchmatrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_hmatrix(pchmatrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_hmatrix_amatrix(pchmatrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_hmatrix(pchmatrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_pre_matrix for sparsematrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_sparsematrix_amatrix(pcsparsematrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_sparsematrix_amatrix(pcsparsematrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_pre_matrix for h2matrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_h2matrix_amatrix(pch2matrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_h2matrix_hmatrix(pch2matrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_h2matrix_amatrix(pch2matrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_h2matrix_hmatrix(pch2matrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_pre_matrix for dh2matrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_dh2matrix_amatrix(pcdh2matrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  LR factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_lr_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_dh2matrix_amatrix(pcdh2matrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$. The matrix @f$B@f$ is given as a
 *  Cholesky factorization and can be applied to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_chol_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_id_pre_matrix for amatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_amatrix(pcamatrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_amatrix_hmatrix(pcamatrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_amatrix(pcamatrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_amatrix_hmatrix(pcamatrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_id_pre_matrix for hmatrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_hmatrix_amatrix(pchmatrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_hmatrix(pchmatrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_hmatrix_amatrix(pchmatrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Hierarchical matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_hmatrix(pchmatrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_id_pre_matrix for sparsematrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_sparsematrix_amatrix(pcsparsematrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_sparsematrix_amatrix(pcsparsematrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A Sparse matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_sparsematrix_hmatrix(pcsparsematrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_id_pre_matrix for h2matrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_h2matrix_amatrix(pch2matrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_h2matrix_hmatrix(pch2matrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_h2matrix_amatrix(pch2matrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_h2matrix_hmatrix(pch2matrix A, pchmatrix chol);

/****************************************************
 * Norm2diff_id_pre_matrix for dh2matrix
 ****************************************************/

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Dense LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_dh2matrix_amatrix(pcdh2matrix A, pcamatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a LR factorization and can be applied to some
 *  vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a LR decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param LR Hierarchical LR factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_lr_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix LR);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Dense Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_dh2matrix_amatrix(pcdh2matrix A, pcamatrix chol);

/** @brief Approximate the spectral norm @f$\|I - B^{-1}A\|_2@f$.
 *  The matrix @f$B@f$ is given as a Cholesky factorization and can be applied
 *  to some vector.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(I - B^{-1}A)^* (I - B^{-1}A)@f$ and computing
 *  the square root of the resulting eigenvalue approximation.
 *  @f$B@f$ is given as a Cholesky decomposition of @f$A@f$.
 *
 *  @param A D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param chol Hierarchical Cholesky factorization of the matrix @f$B@f$.
 *  @returns Approximation of @f$\|I - B^{-1}A\|_2@f$. */
HEADER_PREFIX real
norm2diff_id_chol_dh2matrix_hmatrix(pcdh2matrix A, pchmatrix chol);

/**
 *  @}
 *  */

#endif /* MATRIXNORMS_H_ */
