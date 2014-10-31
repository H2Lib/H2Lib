
/**
    @file h2arith.h
    @author Knut Reimer
    @date 20010-2013
    
    This is the @f$\mathcal{H}^2@f$-matrix-arithmetics module of the H2Lib package.
    It uses the @f$\mathcal{H}^2@f$-matrix-structure of \ref h2matrix.h.
    */

/*@todo extend the function prototypes to general cases*/
#ifndef H2ARITH_H
#define H2ARITH_H

#include "clusteroperator.h"
#include "h2matrix.h"
#include "h2compression.h"
#include "h2update.h"

/** @defgroup h2arith h2arith
 *  @brief Algebraic operations with @f$ \mathcal H^2@f$-matrices.
 *
 *  @attention This is an experimental module, not all functions
 *    have been thoroughly tested, some special cases are not yet
 *    implemented.
 *  @{ */

/* ------------------------------------------------------------
   Auxiliary functions
   ------------------------------------------------------------ */

/** @brief builds a lower triangular H2-matrix based on block tree b
    @return new h2matrix a, upper diagonal blocks are set to 
        @f$ \texttt{a->son} == \texttt{a->f} == \texttt{a->u} = 0 @f$
    @param b : gives the block structure
    @param rb : the row cluster basis of the new h2matrix
    @param cb : the column cluster basis of the new h2matrix
  */
HEADER_PREFIX ph2matrix
build_from_block_lower_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb);

/** @brief builds a upper triangular H2-matrix based on block tree b
    @return new h2matrix a, lower diagonal blocks are set to 
        @f$ \texttt{a->son} == \texttt{a->f} == \texttt{a->u} = 0 @f$
    @param b : gives the block structure
    @param rb : the row cluster basis of the new h2matrix
    @param cb : the column cluster basis of the new h2matrix
  */
HEADER_PREFIX ph2matrix
build_from_block_upper_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb);


/** @brief Clears all upper diagonal blocks of a @ref h2matrix.
 * 
 *  @param a h2matrix.
 */
HEADER_PREFIX void
tolower_h2matrix(ph2matrix a);

/* ------------------------------------------------------------
   Sum
   ------------------------------------------------------------ */

/** @brief @f$ B \gets  B + a @f$
    @param a : constant hmatrix
    @param B : overwritten by @f$ B + a @f$
    @param rwf : has to be the father of the weights of the row clusterbasis of B,\n
        e.g. initialized by prepare_row_clusteroperator
    @param cwf : has to be the father of the weights of the col clusterbasis of B,\n
        e.g. initialized by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    */
HEADER_PREFIX void
add_hmatrix_h2matrix(pchmatrix a,
                     ph2matrix B, pclusteroperator rwf, pclusteroperator cwf,
                     ptruncmode tm, real tol);

/* ------------------------------------------------------------
   Convert functions
   ------------------------------------------------------------ */

/** @brief converts an h2matrix to a new amatrix
    @return new amatrix @f$ B = a @f$
    @param atrans : set if a has to be transposed
    @param a : constant h2matrix
    @todo replace by add_h2matrix_amatrix
    */
HEADER_PREFIX pamatrix
convert_h2matrix_amatrix(bool atrans, pch2matrix a);

/* static or in clusterbasis ??? */
/** @brief computes
    @f$ Y = cb\rightarrow V \cdot Yt_{|cb\rightarrow k \times Yt\rightarrow cols} + Y @f$
    */
HEADER_PREFIX void
fastaddmul_clusterbasis_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Y);

/** @brief converts an uniform matrix @f$ a @f$ to a new rkmatrix \n
    uses fastaddmul_clusterbasis_amatrix
    @return new rkmatrix @f$ B = a @f$
    @param atrans set if @f$ a @f$ has to be transposed, i.e. @f$ B = a^* @f$
    @param a constant h2matrix
    @todo replace by add_uniform_rkmatrix
    */
HEADER_PREFIX prkmatrix
convert_uniform_rkmatrix(bool atrans, pcuniform a);

/** @brief converts an h2matrix to a new hmatrix \n
    @return new hmatrix @f$ B = a @f$, \n
        B has the same block structure as a
    @param a constant h2matrix
    @todo replace by add_h2matrix_hmatrix
    */
HEADER_PREFIX phmatrix
convert_h2matrix_hmatrix(pch2matrix a);

/* ------------------------------------------------------------
   Multiplication
   ------------------------------------------------------------ */

/** @brief computes @f$ a \cdot B @f$ and converts it to a new rkmatrix
    @return new rkmatrix @f$ C = a \cdot B @f$
    @param a : constant h2matrix
    @param btrans : set if B has to be transposed
    @param B : constant h2matrix
    @param tol : tolerance of truncation
    */
HEADER_PREFIX prkmatrix
mul_h2matrix_rkmatrix(pch2matrix a, bool btrans, pch2matrix B, real tol);

/** @brief @f$ C \gets C + alpha * a * B @f$
    @param alpha : field
    @param a : constant h2matrix
    @param btrans : set if B has to be transposed
    @param B : constant h2matrix
    @param C : overwritten by @f$ C + alpha * a * B @f$
    @param rwf : has to be the father of the weights of the row clusterbasis of C,\n
        e.g. initialised by prepare_row_clusteroperator
    @param cwf : has to be the father of the weights of the col clusterbasis of C,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    */
HEADER_PREFIX void
addmul_h2matrix(field alpha, pch2matrix a, bool btrans, pch2matrix B,
                              ph2matrix C, pclusteroperator rwf, pclusteroperator cwf,
                              ptruncmode tm, real tol);

/** @brief computes @f$ a \cdot B @f$ and converts it to a new h2matrix
    @return new h2matrix @f$ C = alpha \cdot a \cdot B @f$
    @param alpha : field
    @param a : constant h2matrix
    @param B : constant h2matrix
    @param h2 : the new h2matrix has the block structure as h2
    @param tm : options of truncation
    @param tol : tolerance of truncation
    */
HEADER_PREFIX ph2matrix
mul_h2matrix(field alpha, ph2matrix a, ph2matrix B,
             ph2matrix h2, ptruncmode tm, real tol);

/* ------------------------------------------------------------
   Inversion
   ------------------------------------------------------------ */

/** @brief @f$ B \gets a^{-1} @f$
    @param a : original matrix, overwritten by auxilliary results
    @param arwf : has to be the father of the weights of the row clusterbasis of a,\n
        e.g. initialised by prepare_row_clusteroperator
    @param acwf : has to be the father of the weights of the col clusterbasis of a,\n
        e.g. initialised by prepare_col_clusteroperator
    @param B : has to have zero entries and the same block tree structure as a, \n
        e.g. initialised by build_fromcluster_clusterbasis, buildfromblock_h2matrix and clear_h2matrix
    @param brwf : has to be the father of the weights of the row clusterbasis of B,\n
        e.g. initialised by prepare_row_clusteroperator
    @param bcwf : has to be the father of the weights of the col clusterbasis of B,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation 
    */
HEADER_PREFIX void
invert_h2matrix(ph2matrix a, pclusteroperator arwf, pclusteroperator acwf,
                ph2matrix B, pclusteroperator brwf, pclusteroperator bcwf,
                ptruncmode tm, real tol);

/* ------------------------------------------------------------
   Forward / backward substitution
   ------------------------------------------------------------ */

/** @brief computes @f$ x \gets T^{-1} \cdot x @f$, using forward substitution
    @param unit : set if T has unit diagonal
    @param a : constant h2matrix
    @param atrans : set if T is the transposed upper triangular part of a,\n
        else T is the lower triangular part of a
    @param x : overwritten by @f$ T^{-1} \cdot x @f$
    */
HEADER_PREFIX void
lowersolve_h2matrix_avector(bool unit, bool atrans, pch2matrix a, pavector x);

/** @brief computes @f$ x \gets T^{-1} \cdot x @f$, using backward substitution
    @param unit : set if T has unit diagonal
    @param a : constant h2matrix
    @param atrans : set if T is the transposed lower triangular part of a,\n
        else T is the upper triangular part of a
    @param x : overwritten by @f$ T^{-1} \cdot x @f$
    */
HEADER_PREFIX void
uppersolve_h2matrix_avector(bool unit, bool atrans, pch2matrix a, pavector x);



/** @brief computes @f$ X \gets T^{-1} \cdot X @f$ or @f$ X^* \gets T^{-1} \cdot X^* @f$, using forward substitution
    @param unit : set if T has unit diagonal
    @param atrans : set if T is the transposed upper triangular part of a,\n
        else T is the lower triangular part of a
    @param a : constant h2matrix
    @param xtrans : set if X has to be transposed
    @param X : overwritten by @f$ T^{-1} \cdot X @f$ or @f$ X \cdot T^{-*}  @f$
    */
HEADER_PREFIX void
lowersolve_h2matrix_amatrix(bool unit, bool atrans, pch2matrix a,
                            bool xtrans, pamatrix X);


/** @brief computes @f$ X \gets T^{-1} \cdot X @f$ or @f$ X^* \gets T^{-1} \cdot X^* @f$, using backward substitution
    @param unit : set if T has unit diagonal
    @param atrans : set if T is the transposed lower triangular part of a,\n
        else T is the upper triangular part of a
    @param a : constant h2matrix
    @param xtrans : set if X has to be transposed
    @param X : overwritten by @f$ T^{-1} \cdot X @f$
    @todo not tested yet
    */
HEADER_PREFIX void
uppersolve_h2matrix_amatrix(bool unit, bool atrans, pch2matrix a,
                            bool xtrans, pamatrix X);

/** @brief computes @f$ Y \gets T^{-1} \cdot X @f$ or @f$ Y^* \gets T^{-1} \cdot X^* @f$, using forward substitution
    @param aunit : set if T has unit diagonal
    @param atrans : set if T is the transposed upper triangular part of a,\n
        else T is the lower triangular part of a
    @param a : constant amatrix
    @param xytrans : set if X and Y has to be transposed
    @param X : constant h2matrix
    @param Y : overwritten by @f$ T^{-1} \cdot X @f$ or @f$ X \cdot T^{-*}  @f$
    @param rwf : has to be the father of the weights of the row clusterbasis of Y,\n
        e.g. initialised by prepare_row_clusteroperator
    @param cwf : has to be the father of the weights of the col clusterbasis of Y,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    @todo ((atrans == true) && (xytrans == false)) is not implemented yet
    */
HEADER_PREFIX void
lowersolve_amatrix_h2matrix(bool aunit, bool atrans, pcamatrix a,
                            bool xytrans, pch2matrix X,
                            ph2matrix Y, pclusteroperator rwf, pclusteroperator cwf,
                            ptruncmode tm, real tol);

/** @brief computes @f$ Y \gets T^{-1} \cdot Y @f$ or @f$ Y^* \gets T^{-1} \cdot Y^* @f$, using backward substitution
    @param aunit : set if T has unit diagonal
    @param atrans : set if T is the transposed lower triangular part of a,\n
        else T is the upper triangular part of a
    @param a : constant amatrix
    @param ytrans : set if Y has to be transposed
    @param Y : overwritten by @f$ T^{-1} \cdot Y @f$ or @f$ Y \cdot T^{-*}  @f$
    @param rwf : has to be the father of the weights of the row clusterbasis of Y,\n
        e.g. initialised by prepare_row_clusteroperator
    @param cwf : has to be the father of the weights of the col clusterbasis of Y,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    @todo not tested yet
    @todo ((atrans == false) && (xytrans == true)) is not implemented yet
    @todo ((atrans == true) && (xytrans == false)) is not implemented yet
    */
HEADER_PREFIX void
uppersolve_amatrix_h2matrix(bool aunit, bool atrans, pcamatrix a,
                            bool ytrans, ph2matrix Y, pclusteroperator rwf, pclusteroperator cwf,
                            ptruncmode tm, real tol);

/** @brief computes @f$ Y \gets T^{-1} \cdot X @f$ or @f$ Y^T \gets T^{-1} \cdot X^T @f$,
      using forward substitution 
    @param aunit : set if T has unit diagonal
    @param atrans : set if T is the transposed upper triangular part of a,\n
        else T is the lower triangular part of a
    @param a : constant h2matrix
    @param xytrans : set if X and Y has to be transposed
    @param X : overwritten by auxilliary results
    @param xrwf : has to be the father of the weights of the row clusterbasis of X,\n
        e.g. initialised by prepare_row_clusteroperator
    @param xcwf : has to be the father of the weights of the col clusterbasis of X,\n
        e.g. initialised by prepare_col_clusteroperator
    @param Y : overwritten by @f$ T^{-1} \cdot X @f$ or @f$ X \cdot T^{-*} @f$
    @param yrwf : has to be the father of the weights of the row clusterbasis of R,\n
        e.g. initialised by prepare_row_clusteroperator
    @param ycwf : has to be the father of the weights of the col clusterbasis of R,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    @todo ((atrans == true) && (xtrans == false)) is not implemented yet
    */
HEADER_PREFIX void
lowersolve_h2matrix(bool aunit, bool atrans, pch2matrix a,
  bool xytrans, ph2matrix X, pclusteroperator xrwf, pclusteroperator xcwf,
  ph2matrix Y, pclusteroperator yrwf, pclusteroperator ycwf,
  ptruncmode tm, real tol);

/** @brief computes @f$ Y \gets T^{-1} \cdot Y @f$ or @f$ Y^T \gets T^{-1} \cdot Y^T @f$,
      using backward substitution 
    @param aunit : set if T has unit diagonal
    @param atrans : set if T is the transposed lower triangular part of a,\n
        else T is the upper triangular part of a
    @param a : constant h2matrix
    @param ytrans : set if Y has to be transposed
    @param Y : overwritten by @f$ T^{-1} \cdot Y @f$ or @f$ Y \cdot T^{-*} @f$
    @param yrwf : has to be the father of the weights of the row clusterbasis of R,\n
        e.g. initialised by prepare_row_clusteroperator
    @param ycwf : has to be the father of the weights of the col clusterbasis of R,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation
    @todo not tested yet
    @todo ((atrans == true) && (xytrans == false)) is not implemented yet
    */
HEADER_PREFIX void
uppersolve_h2matrix(bool aunit, bool atrans, pch2matrix a,
  bool ytrans, ph2matrix Y, pclusteroperator yrwf, pclusteroperator ycwf,
  ptruncmode tm, real tol);

/* ------------------------------------------------------------
   LR decomposition
   ------------------------------------------------------------ */

/** @brief computes the LR-decomposition @f$ L \cdot R \gets a @f$ of a, L has unit diagonal
    @param a : original matrix, overwritten by auxilliary results
    @param arwf : has to be the father of the weights of the row clusterbasis of a,\n
        e.g. initialised by prepare_row_clusteroperator
    @param acwf : has to be the father of the weights of the column clusterbasis of a,\n
        e.g. initialised by prepare_col_clusteroperator
    @param L : cleared (lower) H2-matrix, has to have the same block tree structure as a,\n
        e.g. initialised by build_fromcluster_clusterbasis and build_from_block_lower_h2matrix
    @param lrwf : has to be the father of the weights of the row clusterbasis of L,\n
        e.g. initialised by prepare_row_clusteroperator
    @param lcwf : has to be the father of the weights of the column clusterbasis of L,\n
        e.g. initialised by prepare_col_clusteroperator
    @param R : cleared (upper) H2-matrix, has to have the same block tree structure as a,\n
        e.g. initialised by build_fromcluster_clusterbasis and build_from_block_upper_h2matrix
    @param rrwf : has to be the father of the weights of the row clusterbasis of L,\n
        e.g. initialised by prepare_row_clusteroperator
    @param rcwf : has to be the father of the weights of the column clusterbasis of R,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation */
HEADER_PREFIX void
lrdecomp_h2matrix(ph2matrix a, pclusteroperator arwf, pclusteroperator acwf,
                  ph2matrix L, pclusteroperator lrwf, pclusteroperator lcwf,
                  ph2matrix R, pclusteroperator rrwf, pclusteroperator rcwf,
                  ptruncmode tm, real tol);

/** @brief @f$ x \gets R^{-1} \cdot L^{-1} \cdot x @f$
    @param L : L is the lower triangular part with unit diagonal
    @param R: R is the upper triangular part.
    @param x : overwritten by @f$ R^{-1} \cdot L^{-1} \cdot x @f$
    */
HEADER_PREFIX void
lrsolve_h2matrix_avector(pch2matrix L, pch2matrix R, pavector x);


/* ------------------------------------------------------------
   Cholesky decomposition
   ------------------------------------------------------------ */

/** @brief initialises a, rwf and cwf for choldecomp_h2matrix
    @param a : upper diagonal blocks are set to 
        @f$ \texttt{a->son} == \texttt{a->f} == \texttt{a->u} = 0 @f$
    @param prwf : *prwf is initialised by prepare_row_clusteroperator for new a
    @param pcwf : *pcwf is initialised by prepare_col_clusteroperator for new a
    @param tm : used in prepare_row/col_clusteroperator
    */
HEADER_PREFIX void
init_cholesky_h2matrix(ph2matrix a, pclusteroperator *prwf, pclusteroperator *pcwf,
									ptruncmode tm);

/** @brief computes the Cholesky decomposition @f$ L \cdot L^* \gets@f$ of a and stores the result in L and R
    @param a : original matrix, overwritten by auxilliary results,\n
        e.g. initialised by init_cholesky_h2matrix
    @param arwf : has to be the father of the weights of the row clusterbasis of a,\n
        e.g. initialised by init_cholesky_h2matrix
    @param acwf : has to be the father of the weights of the column clusterbasis of a,\n
        e.g. initialised by init_cholesky_h2matrix
    @param L : cleared (lower) H2-matrix, has to have the same block tree structure as a,\n
        e.g. initialised by build_fromcluster_clusterbasis and build_from_block_lower_h2matrix
    @param lrwf : has to be the father of the weights of the row clusterbasis of L,\n
        e.g. initialised by prepare_row_clusteroperator
    @param lcwf : has to be the father of the weights of the column clusterbasis of L,\n
        e.g. initialised by prepare_col_clusteroperator
    @param tm : options of truncation
    @param tol : tolerance of truncation */
HEADER_PREFIX void
choldecomp_h2matrix(ph2matrix a, pclusteroperator arwf, pclusteroperator acwf,
                  ph2matrix L, pclusteroperator lrwf, pclusteroperator lcwf,
									ptruncmode tm, real tol);

/** @brief @f$ x \gets L^{-*} \cdot L^{-1} \cdot x @f$
    @param a : L is the lower triangular part of a with non-unit diagonal
    @param x : overwritten by @f$ L^{-*} \cdot L^{-1} \cdot x @f$
    */
HEADER_PREFIX void
cholsolve_h2matrix_avector(pch2matrix a, pavector x);

/** @} */

#endif
