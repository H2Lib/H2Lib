/* ------------------------------------------------------------
 * This is the file "harith.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2014
 * ------------------------------------------------------------ */

/** @file harith.h
 @author Steffen B&ouml;rm
 */

#ifndef HARITH_H
#define HARITH_H

#include "hmatrix.h"
#include "truncation.h"

/** @defgroup harith harith
 *  @brief Algebraic operations with hierarchical matrices.
 *  @{ */

/* ------------------------------------------------------------
 Truncation of an rkmatrix
 ------------------------------------------------------------ */

/** @brief Truncate an rkmatrix,
 *  @f$A \gets \operatorname{trunc}(A,\epsilon)@f$.
 *
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param r Source matrix, will be overwritten by truncated matrix. */
HEADER_PREFIX void
trunc_rkmatrix(pctruncmode tm, real eps, prkmatrix r);

/* ------------------------------------------------------------
 * Matrix addition
 * ------------------------------------------------------------ */

/** @brief Truncated addition of a matrix to a low-rank matrices
 *  in @ref rkmatrix representation,
 *  @f$B \gets \operatorname{trunc}(B+\alpha A,\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix_rkmatrix(field alpha, bool atrans, pcamatrix a, pctruncmode tm,
    real eps, prkmatrix b);

///* @brief Truncated addition of a full matrix to a low-rank matrix
// *  in @ref rkmatrix representation,
// *  @f$B \gets \operatorname{trunc}(B+\alpha A,\epsilon)@f$.
// *
// *  This version uses a rank-revealing QR factorization instead
// *  of the standard singular value decomposition. This method
// *  may lead to larger ranks, but may also be faster in some
// *  applications.
// *
// *  @param alpha Scaling factor @f$\alpha@f$.
// *  @param a Source matrix @f$A@f$.
// *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
// *  @param tm Truncation mode.
// *  @param eps Truncation accuracy @f$\epsilon@f$.
// *  @param b Target matrix @f$B@f$. */
//HEADER_PREFIX void
//add2_amatrix_rkmatrix(field alpha, bool atrans, pcamatrix a,
//		      pctruncmode tm, real eps, prkmatrix b);

/** @brief Addition of a low-rank matrix to a matrix,
 *  @f$B \gets B+\alpha A@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_rkmatrix_amatrix(field alpha, bool atrans, pcrkmatrix a, pamatrix b);

/** @brief Addition of a hierarchical matrix to a matrix,
 *  @f$B \gets B+\alpha A@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_hmatrix_amatrix(field alpha, bool atrans, pchmatrix a, pamatrix b);

/** @brief Truncated addition of two low-rank matrices in @ref rkmatrix
 *  representation,
 *  @f$A \gets \operatorname{trunc}(A+\alpha B,\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param src Source matrix @f$B@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param trg Target matrix @f$A@f$. */
HEADER_PREFIX void
add_rkmatrix(field alpha, pcrkmatrix src, pctruncmode tm, real eps,
    prkmatrix trg);

///* @brief Truncated addition of two low-rank matrices in @ref rkmatrix
// *  representation,
// *  @f$B \gets \operatorname{trunc}(B+\alpha A,\epsilon)@f$.
// *
// *  This version uses a rank-revealing QR factorization instead
// *  of the standard singular value decomposition. This method
// *  may lead to larger ranks, but may also be faster in some
// *  applications.
// *
// *  @param alpha Scaling factor @f$\alpha@f$.
// *  @param a Source matrix @f$A@f$.
// *  @param tm Truncation mode.
// *  @param eps Truncation accuracy @f$\epsilon@f$.
// *  @param b Target matrix @f$B@f$. */
//HEADER_PREFIX void
//add2_rkmatrix(field alpha, pcrkmatrix a,
//	      pctruncmode tm, real eps, prkmatrix b);

/** @brief Locally truncated addition of a matrix in @ref amatrix
 *  representation to a hierarchical matrix in @ref hmatrix representation,
 *  @f$B \gets \operatorname{blocktrunc}(B + \alpha A,\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix_hmatrix(field alpha, bool atrans, pcamatrix a, pctruncmode tm,
    real eps, phmatrix b);

/** @brief Locally truncated addition of a matrix in @ref amatrix
 *  representation to the lower triangular part of a hierarchical matrix
 *  in @ref hmatrix representation,
 *  @f$\operatorname{lower}(B) \gets \operatorname{blocktrunc}(\operatorname{lower}(B + \alpha A),\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_lower_amatrix_hmatrix(field alpha, bool atrans, pcamatrix a, pctruncmode tm,
    real eps, phmatrix b);

/** @brief Locally truncated addition of a low-rank matrix in @ref rkmatrix
 *  representation to a hierarchical matrix in @ref hmatrix representation,
 *  @f$B \gets \operatorname{blocktrunc}(B + \alpha A,\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Low-rank matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_rkmatrix_hmatrix(field alpha, pcrkmatrix a, pctruncmode tm, real eps,
    phmatrix b);

/** @brief Locally truncated addition of a low-rank matrix in @ref rkmatrix
 *  representation to the lower triangular part of a hierarchical matrix
 *  in @ref hmatrix representation,
 *  @f$\operatorname{lower}(B) \gets \operatorname{blocktrunc}(\operatorname{lower}(B + \alpha A),\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Low-rank matrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_lower_rkmatrix_hmatrix(field alpha, pcrkmatrix a, pctruncmode tm, real eps,
    phmatrix b);

/** @brief Locally truncated addition of a hierarchical matrix to
 *  a hierarchical matrix representation,
 *  @f$B \gets \operatorname{blocktrunc}(B + \alpha A,\epsilon)@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a @ref hmatrix @f$A@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param b Target @ref hmatrix @f$B@f$. */
HEADER_PREFIX void
add_hmatrix(field alpha, pchmatrix a, pctruncmode tm, real eps, phmatrix b);

/* ------------------------------------------------------------
 * Splitting and merging H-matrices
 * ------------------------------------------------------------ */

///* @brief Split an @ref amatrix into submatrices.
// *
// *  @param f Original matrix.
// *  @param rc Row cluster.
// *  @param cc Column cluster.
// *  @param rsplit Set to split in row direction.
// *  @param csplit Set to split in column direction.
// *  @param copy Set to copy values into new matrix, otherwise it
// *     will be zero.
// *  @returns Block matrix of depth one. */
//HEADER_PREFIX phmatrix
//split_amatrix(pcamatrix f, pccluster rc, pccluster cc,
//	      bool rsplit, bool csplit, bool copy);
/** @brief Split an @ref amatrix into submatrices referencing the
 *  original matrix.
 *
 *  @param f Original matrix.
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param rsplit Set to split in row direction.
 *  @param csplit Set to split in column direction.
 *  @returns Block matrix of depth one. */
HEADER_PREFIX phmatrix
split_sub_amatrix(pcamatrix f, pccluster rc, pccluster cc, bool rsplit,
    bool csplit);

/** @brief Split an @ref rkmatrix into submatrices.
 *
 *  @param r Original matrix.
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param rsplit Set to split in row direction.
 *  @param csplit Set to split in column direction.
 *  @param copy Set to copy values into new matrix, otherwise it
 *     will be zero.
 *  @returns Block matrix of depth one. */
HEADER_PREFIX phmatrix
split_rkmatrix(pcrkmatrix r, pccluster rc, pccluster cc, bool rsplit,
    bool csplit, bool copy);

/** @brief Split an @ref rkmatrix into submatrices referencing the
 *  original matrix.
 *
 *  @attention Since the blocks of the result contain only references
 *  to the original matrix, they cannot be resized, and changing
 *  them changes the other submatrices in the same row and column.
 *  It's probably best to treat the submatrices as constant.
 *
 *  @param r Original matrix.
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param rsplit Set to split in row direction.
 *  @param csplit Set to split in column direction.
 *  @returns Block matrix of depth one. */
HEADER_PREFIX phmatrix
split_sub_rkmatrix(pcrkmatrix r, pccluster rc, pccluster cc, bool rsplit,
    bool csplit);

///* @brief Split a leaf @ref hmatrix into submatrices.
// *
// *  @param hm Leaf @ref hmatrix.
// *  @param rsplit Set to split in row direction.
// *  @param csplit Set to split in column direction.
// *  @param copy Set to copy values into new matrix, otherwise it
// *     will be zero.
// *  @returns Block matrix of depth one or zero according to
// *     <tt>rsplit</tt> and <tt>csplit</tt>. */
//HEADER_PREFIX phmatrix
//split_hmatrix(phmatrix hm, bool rsplit, bool csplit, bool copy);
//
///* @brief Split a leaf @ref hmatrix into submatrices referencing
// *  the original matrix.
// *
// *  @attention Since the blocks of the result contain only references
// *  to the original matrix, they cannot be resized, and changing
// *  them changes the original matrix.
// *
// *  @param hm Leaf @ref hmatrix.
// *  @param rsplit Set to split in row direction.
// *  @param csplit Set to split in column direction.
// *  @returns Block matrix of depth one or zero according to
// *     <tt>rsplit</tt> and <tt>csplit</tt>. */
//HEADER_PREFIX phmatrix
//split_sub_hmatrix(phmatrix hm,
//		  bool rsplit, bool csplit);

/** @brief Block-merge two rk matrices into a new rk matrix,
 *  @f$A \gets \operatorname{trunc}\left(\begin{pmatrix} A & B \end{pmatrix},\epsilon\right)@f$ or
 *  @f$A \gets \operatorname{trunc}\left(\begin{pmatrix} A\\ B \end{pmatrix},\epsilon\right)@f$.
 *
 *  If <tt>colmerge</tt> is not set, @f$A@f$ and @f$B@f$ have to have the
 *  same number of rows, otherwise the same number of columns.
 *
 *  @attention The matrix @f$A@f$ is the last parameter, although it is
 *  placed in the first row or column.
 *
 *  @param colmerge Set to merge into a column instead of a row.
 *  @param src Source matrix @f$B@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @param trg Target matrix @f$A@f$, will be overwritten by result. */
HEADER_PREFIX void
merge_rkmatrix(bool colmerge, pcrkmatrix src, pctruncmode tm, real eps,
    prkmatrix trg);

/** @brief Merge a depth-one @ref hmatrix containing only
 *  @ref rkmatrix submatrices into a new @ref rkmatrix,
 *  @f$A \gets \operatorname{succtrunc}\left(\begin{pmatrix}
 *  B_{11} & \ldots & B_{1m}\\ \vdots & \ddots & \vdots\\
 *  B_{n1} & \ldots & B_{nm} \end{pmatrix},\epsilon\right)@f$.
 *
 *  @param s Source matrix @f$B@f$
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy @f$\epsilon@f$.
 *  @returns @ref rkmatrix approximating the block matrix. */
HEADER_PREFIX prkmatrix
merge_hmatrix_rkmatrix(pchmatrix s, pctruncmode tm, real eps);

/* ------------------------------------------------------------
 * Matrix multiplication
 * ------------------------------------------------------------ */

/** @brief Multiply a matrix by an @ref rkmatrix,
 *  @f$C \gets C + \alpha A B@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Low-rank matrix @f$A@f$.
 *  @param btrans Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param b Matrix @f$B@f$.
 *  @param ctrans Set if @f$C^*@f$ is to be used instead of @f$C@f$.
 *  @param c Target matrix @f$C@f$. */
HEADER_PREFIX void
addmul_rkmatrix_amatrix_amatrix(field alpha, bool atrans, pcrkmatrix a,
    bool btrans, pcamatrix b, bool ctrans, pamatrix c);

/** @brief Multiply a matrix by a @ref hmatrix,
 *  @f$C \gets C + \alpha A B@f$.
 *
 *  @remark The matrices @f$C@f$ and @f$B@f$ are assumed to be
 *  ordered matching the row and column cluster trees of @f$A@f$,
 *  respectively.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Hierarchical matrix @f$A@f$.
 *  @param btrans Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param bp Matrix @f$B@f$.
 *  @param ctrans Set if @f$C^*@f$ is to be used instead of @f$C@f$.
 *  @param cp Target matrix @f$C@f$. */
HEADER_PREFIX void
addmul_hmatrix_amatrix_amatrix(field alpha, bool atrans, pchmatrix a,
    bool btrans, pcamatrix bp, bool ctrans, pamatrix cp);

/** @brief Multiply two H-matrices,
 *  @f$Z \gets \operatorname{succtrunc}(Z + \alpha X Y,\epsilon)@f$.
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
addmul_hmatrix(field alpha, bool xtrans, pchmatrix x, bool ytrans, pchmatrix y,
    pctruncmode tm, real eps, phmatrix z);

/** @brief Multiply two H-matrices, computing only the lower triangular
 *  part of the result,
 *  @f$\operatorname{lower}(Z) \gets \operatorname{succtrunc}(\operatorname{lower}(Z + \alpha X Y),\epsilon)@f$.
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
addmul_lower_hmatrix(field alpha, bool xtrans, pchmatrix x, bool ytrans,
    pchmatrix y, pctruncmode tm, real eps, phmatrix z);

/* ------------------------------------------------------------
 * Matrix inversion
 * ------------------------------------------------------------ */

/** @brief Invert an H-matrix,
 *  @f$A \gets \operatorname{succtrunc}(A^{-1},\epsilon)@f$.
 *
 *  @param a Matrix @f$A@f$, gets overwritten with @f$A^{-1}@f$.
 *  @param work Auxiliary storage, has to be zero and be of the
 *     same structure as <tt>a</tt>.
 *     Consider using @ref clonestructure_hmatrix.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
invert_hmatrix(phmatrix a, phmatrix work, pctruncmode tm, real eps);

/* ------------------------------------------------------------
 * Triangular matrices
 * ------------------------------------------------------------ */

/** @brief Create a duplicate of the lower triangular part of an
 *  H-matrix.
 *
 *  @param aunit Set if the duplicate has to have unit diagonal.
 *  @param a Source matrix @f$A@f$.
 *  @returns Lower triangular matrix containing the lower triangular
 *     part of @f$A@f$, with unit diagonal if <tt>aunit==true</tt>. */
HEADER_PREFIX phmatrix
clone_lower_hmatrix(bool aunit, pchmatrix a);

/** @brief Create a duplicate of the upper triangular part of an
 *  H-matrix.
 *
 *  @param aunit Set if the duplicate has to have unit diagonal.
 *  @param a Source matrix @f$A@f$.
 *  @returns Upper triangular matrix containing the upper triangular
 *     part of @f$A@f$, with unit diagonal if <tt>aunit==true</tt>. */
HEADER_PREFIX phmatrix
clone_upper_hmatrix(bool aunit, pchmatrix a);

/* ------------------------------------------------------------
 * Solve triangular systems
 * ------------------------------------------------------------ */

/** @brief Solve a triangular system,
 *  @f$x \gets R^{-1} x@f$ or @f$x \gets R^{-*} x@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to be upper triangular.
 *  @param aunit Set if @f$R@f$ has unit diagonal.
 *  @param atrans Set if @f$R^*@f$ instead of @f$R@f$ is to be used.
 *  @param a Matrix containing the upper triangular part of @f$R@f$.
 *     The strictly lower triangular part is not used.
 *  @param xp Input vector @f$x@f$, will be overwritten by the solution.
 *     The vector is assumed to be in cluster numbering according
 *     to <tt>a->cc</tt> respectively <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularinvmul_hmatrix_avector(bool alower, bool aunit, bool atrans,
    pchmatrix a, pavector xp);

/** @brief Solve a triangular system @f$x \gets A^{-1} x@f$ or
 *  @f$x \gets (A^*)^{-1} x@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to be upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param x Input vector @f$x@f$, will be overwritten by the solution. */
HEADER_PREFIX void
triangularsolve_hmatrix_avector(bool alower, bool aunit, bool atrans,
    pchmatrix a, pavector x);

/** @brief Solve a triangular system @f$X \gets A^{-1} X@f$ or
 *  @f$X \gets (A^*)^{-1} X@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to by upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xp Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularinvmul_hmatrix_amatrix(bool alower, bool aunit, bool atrans,
    pchmatrix a, bool xtrans, pamatrix xp);

/** @brief Solve a triangular system,
 *  @f$X \gets \operatorname{succtrunc}(A^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(A^{-*} X, \epsilon)@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to be upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the triangular part of @f$A@f$.
 *     The remainder is not used.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xp Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularinvmul_hmatrix(bool alower, bool aunit, bool atrans, pchmatrix a,
    pctruncmode tm, real eps, bool xtrans, phmatrix xp);

/** @brief Solve a triangular system,
 *  @f$X \gets \operatorname{succtrunc}(A^{-1} X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(A^{-*} X, \epsilon)@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to by upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xp Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularinvmul_amatrix_hmatrix(bool alower, bool aunit, bool atrans,
    pcamatrix a, pctruncmode tm, real eps, bool xtrans, phmatrix xp);

/* ------------------------------------------------------------
 * Evaluate triangular matrices
 * ------------------------------------------------------------ */

/** @brief Evaluate a triangular system @f$x \gets A x@f$ or
 *  @f$x \gets A^* x@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to be upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param xp Input vector @f$x@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularmul_hmatrix_avector(bool alower, bool aunit, bool atrans, pchmatrix a,
    pavector xp);

/** @brief Evaluate a triangular system @f$x \gets A x@f$ or
 *  @f$x \gets A^* x@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to be upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param x Input vector @f$x@f$, will be overwritten by the solution. */
HEADER_PREFIX void
triangulareval_hmatrix_avector(bool alower, bool aunit, bool atrans,
    pchmatrix a, pavector x);

/** @brief Evaluate a triangular system @f$X \gets A X@f$ or
 *  @f$X \gets A^* X@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to by upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xp Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularmul_hmatrix_amatrix(bool alower, bool aunit, bool atrans, pchmatrix a,
    bool xtrans, pamatrix xp);

/** @brief Evaluate a triangular system,
 *  @f$X \gets \operatorname{succtrunc}(A X, \epsilon)@f$ or
 *  @f$X \gets \operatorname{succtrunc}(A^{*} X, \epsilon)@f$.
 *
 *  @param alower Set if @f$A@f$ is lower triangular, otherwise it
 *     is assumed to by upper triangular.
 *  @param aunit Set if @f$A@f$ has unit diagonal.
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the lower or upper triangular part of @f$A@f$.
 *     The remainder of the matrix is not used.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @param xtrans set if @f$X^*@f$ instead of @f$X@f$ is to be used.
 *  @param xp Input matrix @f$X@f$, will be overwritten by the solution.
 *     The rows are assumed to be in cluster numbering according
 *     to <tt>a->rc</tt>. */
HEADER_PREFIX void
triangularmul_hmatrix(bool alower, bool aunit, bool atrans, pchmatrix a,
    pctruncmode tm, real eps, bool xtrans, phmatrix xp);

/* ------------------------------------------------------------
 * Triangular factorizations
 * ------------------------------------------------------------ */

/** @brief Compute the LR factorization,
 *  @f$A \approx L R@f$.
 *
 *  The lower triangular part @f$L@f$ has unit diagonal, only its
 *  strict lower triangular part is stored in the strict lower
 *  triangular part of the source matrix.
 *
 *  The upper triangular part is stored in the upper triangular part
 *  of the source matrix.
 *
 *  @param a Source matrix @f$A@f$, will be overwritten
 *     by @f$L@f$ and @f$R@f$.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy. */
HEADER_PREFIX void
lrdecomp_hmatrix(phmatrix a, pctruncmode tm, real eps);

/** @brief Solve the linear systems @f$A x = b@f$
 *  using the LR factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param a Matrix containing the LR factorization in the form returned
 *     by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$b@f$, will be overwritten by
 *     solution vector @f$x@f$. */
HEADER_PREFIX void
lrsolve_n_hmatrix_avector(pchmatrix a, pavector x);

/** @brief Solve the linear systems @f$A^* x = b@f$
 *  using the LR factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param a Matrix containing the LR factorization in the form returned
 *     by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$b@f$, will be overwritten by
 *     solution vector @f$x@f$. */
HEADER_PREFIX void
lrsolve_t_hmatrix_avector(pchmatrix a, pavector x);

/** @brief Solve the linear systems @f$A x = b@f$ or @f$A^* x = b@f$
 *  using the LR factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the LR factorization in the form returned
 *     by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$b@f$, will be overwritten by
 *     solution vector @f$x@f$. */
HEADER_PREFIX void
lrsolve_hmatrix_avector(bool atrans, pchmatrix a, pavector x);

/**
 * @brief Evaluate the linear system @f$x \gets A x@f$ using the LR
 *  factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param a Matrix containing the LR factorization in the form
 *    returned by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$x@f$, that gets overwritten by
 *     the result. */
HEADER_PREFIX void
lreval_n_hmatrix_avector(pchmatrix a, pavector x);

/**
 * @brief Evaluate the linear system @f$x \gets A^* x@f$ using the LR
 *  factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param a Matrix containing the LR factorization in the form
 *    returned by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$x@f$, that gets overwritten by
 *     the result. */
HEADER_PREFIX void
lreval_t_hmatrix_avector(pchmatrix a, pavector x);

/**
 * @brief Evaluate the linear system @f$x \gets A x@f$ or @f$x \gets A^* x@f$
 *  using the LR factorization provided by @ref lrdecomp_hmatrix.
 *
 *  @param atrans Set if @f$A^*@f$ instead of @f$A@f$ is to be used.
 *  @param a Matrix containing the LR factorization in the form
 *    returned by @ref lrdecomp_hmatrix.
 *  @param x Right-hand side vector @f$x@f$, that gets overwritten by
 *     the result. */
HEADER_PREFIX void
lreval_hmatrix_avector(bool atrans, pchmatrix a, pavector x);

/** @brief Compute the Cholesky factorization,
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
choldecomp_hmatrix(phmatrix a, pctruncmode tm, real eps);

/** @brief Solve the linear system @f$A x = b@f$ using the Cholesky
 *  factorization provided by @ref choldecomp_hmatrix.
 *
 *  @param a Matrix containing the Cholesky factorization in the form
 *    returned by @ref choldecomp_hmatrix.
 *  @param x Right-hand side vector @f$b@f$, will be overwritten by
 *     solution vector @f$x@f$. */
HEADER_PREFIX void
cholsolve_hmatrix_avector(pchmatrix a, pavector x);

/**
 * @brief Evaluate the linear system @f$x \gets A x@f$ using the Cholesky
 *  factorization provided by @ref choldecomp_hmatrix.
 *
 *  @param a Matrix containing the Cholesky factorization in the form
 *    returned by @ref choldecomp_hmatrix.
 *  @param x Right-hand side vector @f$x@f$, that gets overwritten by
 *     the result. */
HEADER_PREFIX void
choleval_hmatrix_avector(pchmatrix a, pavector x);

/** @} */

#endif
