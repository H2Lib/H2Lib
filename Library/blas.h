/* ------------------------------------------------------------
 This is the file "blas.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file blas.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef BLAS_H_
#define BLAS_H_

#ifdef USE_BLAS

/* C STD LIBRARY */
/* CORE 0 */
#include "settings.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

/** \defgroup blas blas
 *  @brief This module servers as a wrapper for all calls to BLAS and LAPACK
 *  routines
 *
 *  For more detailed information about the imported functions, please refer
 *  e.g. to http://www.netlib.org/lapack/explore-html/index.html
 *  @{ */

/****************************************************
 * BLAS level 1
 ****************************************************/

/****************************************************
 * _DOT / _RDOT
 ****************************************************/

/**
 * @brief Compute the dot product of two field vectors @f$\vec x@f$
 * and @f$\vec y@f$ of length @p n.
 *
 *   @param n Length of both vectors @p x and @p y.
 *   @param x First field vector of length @p n.
 *   @param incx Stride for elements of @p x.
 *   @param y Second field vector of length @p n.
 *   @param incy Stride for elements of @p y.
 *   @return Dot product of @p x and @p y: @f$ \langle \vec x, \vec y \rangle_2
 *   = \sum_{i=1}^n \bar x_i y_i@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX float _Complex
cdotc_(const unsigned *n, const float _Complex * x, const unsigned *incx,
    const float _Complex * y, const unsigned *incy);
/** @endcond */
#define h2_dot(n, x, incx, y, incy) cdotc_(n, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX float
sdot_(const unsigned *n,
    const float *x, const unsigned *incx, const float *y,
    const unsigned *incy);
/** @endcond */
#define h2_dot(n, x, incx, y, incy) sdot_(n, x, incx, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX double _Complex
zdotc_(const unsigned *n, const double _Complex * x, const unsigned *incx,
    const double _Complex * y, const unsigned *incy);
/** @endcond */
#define h2_dot(n, x, incx, y, incy) zdotc_(n, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX double
ddot_(const unsigned *n, const double *x, const unsigned *incx, const double *y,
    const unsigned *incy);
/** @endcond */
#define h2_dot(n, x, incx, y, incy) ddot_(n, x, incx, y, incy)
#endif
#endif

/**
 * @brief Compute the dot product of two real vectors @f$\vec x@f$
 * and @f$\vec y@f$ of length @p n.
 *
 *   @param n Length of both vectors @p x and @p y.
 *   @param x First real vector of length @p n.
 *   @param incx Stride for elements of @p x.
 *   @param y Second real vector of length @p n.
 *   @param incy Stride for elements of @p y.
 *   @return Dot product of @p x and @p y: @f$ \langle \vec x, \vec y \rangle_2
 *   = \sum_{i=1}^n x_i y_i@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX float
sdot_(const unsigned *n, const float *x, const unsigned *incx, const float *y,
    const unsigned *incy);
/** @endcond */
#define h2_rdot(n, x, incx, y, incy) sdot_(n, x, incx, y, incy)
#else
#define h2_rdot(n, x, incx, y, incy) sdot_(n, x, incx, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX double
ddot_(const unsigned *n, const double *x, const unsigned *incx, const double *y,
    const unsigned *incy);
/** @endcond */
#define h2_rdot(n, x, incx, y, incy) ddot_(n, x, incx, y, incy)
#else
#define h2_rdot(n, x, incx, y, incy) ddot_(n, x, incx, y, incy)
#endif
#endif

/****************************************************
 * _AXPY /_RAXPY
 ****************************************************/

/**
 * @brief Add a field vector @f$\alpha \vec x@f$ to another field
 * vector @f$\vec y@f$.
 *
 * @param n Length of both vectors @p x and @p y.
 * @param alpha Field constant.
 * @param x First field vector of length @p n.
 * @param incx Stride for elements of @p x.
 * @param y Second field vector of length @p n.
 * @param incy Stride for elements of @p y.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
caxpy_(const unsigned *n, const float _Complex * alpha,
    const float _Complex * x, const unsigned *incx, float _Complex * y,
    const unsigned *incy);
/** @endcond */
#define h2_axpy(n, alpha, x, incx, y, incy) caxpy_(n, alpha, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
saxpy_(const unsigned *n, const float * alpha,
    const float * x, const unsigned *incx, float * y, const unsigned *incy);
/** @endcond */
#define h2_axpy(n, alpha, x, incx, y, incy) saxpy_(n, alpha, x, incx, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zaxpy_(const unsigned *n, const double _Complex *alpha,
    const double _Complex *x, const unsigned *incx, double _Complex *y,
    const unsigned *incy);
/** @endcond */
#define h2_axpy(n, alpha, x, incx, y, incy) zaxpy_(n, alpha, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
daxpy_(const unsigned *n, const double *alpha, const double *x,
    const unsigned *incx, double *y, const unsigned *incy);
/** @endcond */
#define h2_axpy(n, alpha, x, incx, y, incy) daxpy_(n, alpha, x, incx, y, incy)
#endif
#endif

/**
 * @brief Add a real vector @f$\alpha \vec x@f$ to another field
 * vector @f$\vec y@f$.
 *
 * @param n Length of both vectors @p x and @p y.
 * @param alpha Real constant.
 * @param x First real vector of length @p n.
 * @param incx Stride for elements of @p x.
 * @param y Second real vector of length @p n.
 * @param incy Stride for elements of @p y.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
saxpy_(const unsigned *n, const float * alpha, const float * x,
    const unsigned *incx, float * y, const unsigned *incy);
/** @endcond */
#define h2_raxpy(n, alpha, x, incx, y, incy) saxpy_(n, alpha, x, incx, y, incy)
#else
#define h2_raxpy(n, alpha, x, incx, y, incy) saxpy_(n, alpha, x, incx, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
daxpy_(const unsigned *n, const double *alpha, const double *x,
    const unsigned *incx, double *y, const unsigned *incy);
/** @endcond */
#define h2_raxpy(n, alpha, x, incx, y, incy) daxpy_(n, alpha, x, incx, y, incy)
#else
#define h2_raxpy(n, alpha, x, incx, y, incy) daxpy_(n, alpha, x, incx, y, incy)
#endif
#endif

/****************************************************
 * _SCAL / _RSCAL
 ****************************************************/

/**
 * @brief Scales a field vector @f$\vec x@f$ by a field scalar @f$\alpha@f$.
 *
 * @param n Length of vector @p x.
 * @param alpha Field constant.
 * @param x Field vector of length @p n.
 * @param incx Stride for elements of @p x.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cscal_(const unsigned *n, const float _Complex *alpha, float _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_scal(n, alpha, x, incx) cscal_(n, alpha, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sscal_(const unsigned *n, const float *alpha, float *x, const unsigned *incx);
/** @endcond */
#define h2_scal(n, alpha, x, incx) sscal_(n, alpha, x, incx)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zscal_(const unsigned *n, const double _Complex *alpha, double _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_scal(n, alpha, x, incx) zscal_(n, alpha, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dscal_(const unsigned *n, const double *alpha, double *x, const unsigned *incx);
/** @endcond */
#define h2_scal(n, alpha, x, incx) dscal_(n, alpha, x, incx)
#endif
#endif

/**
 * @brief Scales a field vector @f$\vec x@f$ by a real scalar @f$\alpha@f$.
 *
 * @param n Length of vector @p x.
 * @param alpha Real constant.
 * @param x Field vector of length @p n.
 * @param incx Stride for elements of @p x.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
csscal_(const unsigned *n, const float *alpha, float _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_rscal(n, alpha, x, incx) csscal_(n, alpha, x, incx)
#else
#define h2_rscal(n, alpha, x, incx) sscal_(n, alpha, x, incx)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zdscal_(const unsigned *n, const double *alpha, double _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_rscal(n, alpha, x, incx) zdscal_(n, alpha, x, incx)
#else
#define h2_rscal(n, alpha, x, incx) dscal_(n, alpha, x, incx)
#endif
#endif

/****************************************************
 * _NRM2
 ****************************************************/

/**
 * @brief Computes the 2-norm @f$\|\vec x\|_2@f$ of a field vector @f$\vec x@f$.
 *
 * @param n Length of vector @p x.
 * @param x Field vector of length @p n.
 * @param incx Stride for elements of @p x.
 * @return @f$\|x\|_2@f$.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX float
scnrm2_(const unsigned *n, const float _Complex *x, const unsigned *incx);
/** @endcond */
#define h2_nrm2(n, x, incx) scnrm2_(n, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX float
snrm2_(const unsigned *n, const float *x, const unsigned *incx);
/** @endcond */
#define h2_nrm2(n, x, incx) snrm2_(n, x, incx)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX double
dznrm2_(const unsigned *n, const double _Complex *x, const unsigned *incx);
/** @endcond */
#define h2_nrm2(n, x, incx) dznrm2_(n, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX double
dnrm2_(const unsigned *n, const double *x, const unsigned *incx);
/** @endcond */
#define h2_nrm2(n, x, incx) dnrm2_(n, x, incx)
#endif
#endif

/****************************************************
 * _SWAP
 ****************************************************/

/**
 * @brief Swaps to entries of two field vectors @f$\vec x@f$ and @f$\vec y@f$.
 *
 * @param n Length of both vectors @p x and @p y.
 * @param x First field vector of length @p n.
 * @param incx Stride for elements of @p x.
 * @param y Second field vector of length @p n.
 * @param incy Stride for elements of @p y.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cswap_(const unsigned *n, float _Complex *x, const unsigned *incx,
    float _Complex *y, const unsigned *incy);
/** @endcond */
#define h2_swap(n, x, incx, y, incy) cswap_(n, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sswap_(const unsigned *n, float *x, const unsigned *incx, float *y,
    const unsigned *incy);
/** @endcond */
#define h2_swap(n, x, incx, y, incy) sswap_(n, x, incx, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zswap_(const unsigned *n, double _Complex *x, const unsigned *incx,
    double _Complex *y, const unsigned *incy);
/** @endcond */
#define h2_swap(n, x, incx, y, incy) zswap_(n, x, incx, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dswap_(const unsigned *n, double *x, const unsigned *incx, double *y,
    const unsigned *incy);
/** @endcond */
#define h2_swap(n, x, incx, y, incy) dswap_(n, x, incx, y, incy)
#endif
#endif

/****************************************************
 * BLAS level 2
 ****************************************************/

/****************************************************
 * _GEMV
 ****************************************************/

/**
 * @brief Compute matrix vector product of a general matrix @f$A@f$ with input
 * vector @f$\vec x@f$ and output vector @f$y@f$.
 *
 * The operation @f$\vec y \gets \beta \vec y + \alpha A \vec x@f$ or
 *                @f$\vec y \gets \beta \vec y + \alpha A^T \vec x@f$ or
 *                @f$\vec y \gets \beta \vec y + \alpha \bar A^T \vec x@f$
 * will be performed depending on the value of @p trans and the current type of
 * @ref field.
 *
 * @param trans If set to "Conjugate transposed" @f$\bar A^T@f$ will be used.<br>
 *              If set to "Transposed" @f$A^T@f$ will be used.<br>
 *              In all other cases @f$A@f$ will be used for the computation.
 * @param m Number of rows for the matrix @p a.
 * @param n Number of columns for the matrix @p a.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param x Input field vector @f$\vec x@f$. Dimension of @p x should be at least
 *        @p n.
 * @param incx Stride for elements of @p x.
 * @param beta Field scalar value used in the above mentioned computation.
 * @param y Output field vector @f$\vec y@f$. Dimension of @p y should be at least
 *        @p m.
 * @param incy Stride for elements of @p y.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgemv_(const char *trans, const unsigned *m, const unsigned *n,
    const float _Complex *alpha, const float _Complex *a, const unsigned *lda,
    const float _Complex *x, const unsigned *incx, const float _Complex *beta,
    float _Complex *y, const unsigned *incy);
/** @endcond */
#define h2_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) cgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sgemv_(const char *trans, const unsigned *m, const unsigned *n,
    const float *alpha, const float *a, const unsigned *lda, const float *x,
    const unsigned *incx, const float *beta, float *y, const unsigned *incy);
/** @endcond */
#define h2_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) sgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgemv_(const char *trans, const unsigned *m, const unsigned *n,
    const double _Complex *alpha, const double _Complex *a, const unsigned *lda,
    const double _Complex *x, const unsigned *incx, const double _Complex *beta,
    double _Complex *y, const unsigned *incy);
/** @endcond */
#define h2_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dgemv_(const char *trans, const unsigned *m, const unsigned *n,
    const double *alpha, const double *a, const unsigned *lda, const double *x,
    const unsigned *incx, const double *beta, double *y, const unsigned *incy);
/** @endcond */
#define h2_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
#endif
#endif

/****************************************************
 * _TRMV
 ****************************************************/

/**
 * @brief Compute matrix vector product of a triangular matrix @f$A@f$ with input
 * and output vector @f$\vec x@f$.
 *
 * The operation @f$\vec x \gets A \vec x@f$ or
 *                @f$\vec x \gets A^T \vec x@f$ or
 *                @f$\vec x \gets \bar A^T \vec x@f$
 * will be performed depending on the value of @p trans and the current type of
 * @ref field.
 *
 * @param uplo Determines if an "upper" right or a "lower" left triangular matrix
 *   @f$A@f$ should be used.
 * @param trans If set to "Conjugate transposed" @f$\bar A^T@f$ will be used.<br>
 *              If set to "Transposed" @f$A^T@f$ will be used.<br>
 *              In all other cases @f$A@f$ will be used for the computation.
 * @param diag Determines if a "Unit" diagonal should be used or "Non-unit"
 *        matrix coefficients stored within the matrix @f$A@f$.
 * @param n Number of row and columns for the matrix @p a.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param x Input/output field vector @f$\vec x@f$. Dimension of @p x should be at least
 *        @p n.
 * @param incx Stride for elements of @p x.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ctrmv_(const char *uplo, const char *trans, const char *diag, const unsigned *n,
    const float _Complex *a, const unsigned *lda, float _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_trmv(uplo, trans, diag, n, a, lda, x, incx) ctrmv_(uplo, trans, diag, n, a, lda, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
strmv_(const char *uplo, const char *trans, const char *diag, const unsigned *n,
    const float *a, const unsigned *lda, float *x, const unsigned *incx);
/** @endcond */
#define h2_trmv(uplo, trans, diag, n, a, lda, x, incx) strmv_(uplo, trans, diag, n, a, lda, x, incx)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ztrmv_(const char *uplo, const char *trans, const char *diag, const unsigned *n,
    const double _Complex *a, const unsigned *lda, double _Complex *x,
    const unsigned *incx);
/** @endcond */
#define h2_trmv(uplo, trans, diag, n, a, lda, x, incx) ztrmv_(uplo, trans, diag, n, a, lda, x, incx)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dtrmv_(const char *uplo, const char *trans, const char *diag, const unsigned *n,
    const double *a, const unsigned *lda, double *x, const unsigned *incx);
/** @endcond */
#define h2_trmv(uplo, trans, diag, n, a, lda, x, incx) dtrmv_(uplo, trans, diag, n, a, lda, x, incx)
#endif
#endif

/****************************************************
 * _GER
 ****************************************************/

/**
 * @brief Adds a rank-1-update @f$\vec x \vec{\bar  y}^T@f$ to a matrix @f$A@f$.
 *
 * The Computation @f$A \gets A + \alpha \vec x \vec{\bar  y}^T@f$ is performed.
 *
 * @param m Number of rows for the matrix @p a.
 * @param n Number of columns for the matrix @p a.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param x First field vector @f$\vec x@f$. Dimension of @p x should be @p m.
 * @param incx Stride for elements of @p x.
 * @param y Second field vector @f$\vec y@f$. Dimension of @p y should be @p n.
 * @param incy Stride for elements of @p y.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgerc_(const unsigned *m, const unsigned *n, const float _Complex *alpha,
    const float _Complex *x, const unsigned *incx, const float _Complex *y,
    const unsigned *incy, float _Complex *a, const unsigned *lda);
/** @endcond */
#define h2_ger(m, n, alpha, x, incx, y, incy, a, lda) cgerc_(m, n, alpha, x, incx, y, incy, a, lda)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sger_(const unsigned *m, const unsigned *n, const float *alpha,
    const float *x, const unsigned *incx, const float *y,
    const unsigned *incy, float *a, const unsigned *lda);
/** @endcond */
#define h2_ger(m, n, alpha, x, incx, y, incy, a, lda) sger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgerc_(const unsigned *m, const unsigned *n, const double _Complex *alpha,
    const double _Complex *x, const unsigned *incx, const double _Complex *y,
    const unsigned *incy, double _Complex *a, const unsigned *lda);
/** @endcond */
#define h2_ger(m, n, alpha, x, incx, y, incy, a, lda) zgerc_(m, n, alpha, x, incx, y, incy, a, lda)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dger_(const unsigned *m, const unsigned *n, const double *alpha,
    const double *x, const unsigned *incx, const double *y,
    const unsigned *incy, double *a, const unsigned *lda);
/** @endcond */
#define h2_ger(m, n, alpha, x, incx, y, incy, a, lda) dger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#endif

/**
 * @brief Adds a rank-1-update @f$\vec x \vec{\bar  y}^T@f$ to a matrix @f$A@f$.
 *
 * The Computation @f$A \gets A + \alpha \vec x \vec{\bar  y}^T@f$ is performed.
 *
 * @param m Number of rows for the matrix @p a.
 * @param n Number of columns for the matrix @p a.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param x First field vector @f$\vec x@f$. Dimension of @p x should be @p m.
 * @param incx Stride for elements of @p x.
 * @param y Second field vector @f$\vec y@f$. Dimension of @p y should be @p n.
 * @param incy Stride for elements of @p y.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
#define h2_gerc(m, n, alpha, x, incx, y, incy, a, lda) cgerc_(m, n, alpha, x, incx, y, incy, a, lda)
#else
#define h2_gerc(m, n, alpha, x, incx, y, incy, a, lda) sger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#else
#ifdef USE_COMPLEX
#define h2_gerc(m, n, alpha, x, incx, y, incy, a, lda) zgerc_(m, n, alpha, x, incx, y, incy, a, lda)
#else
#define h2_gerc(m, n, alpha, x, incx, y, incy, a, lda) dger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#endif

/**
 * @brief Adds a rank-1-update @f$\vec x \vec y^T@f$ to a matrix @f$A@f$.
 *
 * The Computation @f$A \gets A + \alpha \vec x \vec y^T@f$ is performed.
 *
 * @param m Number of rows for the matrix @p a.
 * @param n Number of columns for the matrix @p a.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param x First field vector @f$\vec x@f$. Dimension of @p x should be @p m.
 * @param incx Stride for elements of @p x.
 * @param y Second field vector @f$\vec y@f$. Dimension of @p y should be @p n.
 * @param incy Stride for elements of @p y.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgeru_(const unsigned *m, const unsigned *n, const float _Complex *alpha,
    const float _Complex *x, const unsigned *incx, const float _Complex *y,
    const unsigned *incy, float _Complex *a, const unsigned *lda);
/** @endcond */
#define h2_geru(m, n, alpha, x, incx, y, incy, a, lda) cgeru_(m, n, alpha, x, incx, y, incy, a, lda)
#else
#define h2_geru(m, n, alpha, x, incx, y, incy, a, lda) sger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgeru_(const unsigned *m, const unsigned *n, const double _Complex *alpha,
    const double _Complex *x, const unsigned *incx, const double _Complex *y,
    const unsigned *incy, double _Complex *a, const unsigned *lda);
/** @endcond */
#define h2_geru(m, n, alpha, x, incx, y, incy, a, lda) zgeru_(m, n, alpha, x, incx, y, incy, a, lda)
#else
#define h2_geru(m, n, alpha, x, incx, y, incy, a, lda) dger_(m, n, alpha, x, incx, y, incy, a, lda)
#endif
#endif

/****************************************************
 * _SYR
 ****************************************************/

/**
 * @brief Adds a symmetric/hermetian rank-1-update
 * @f$\vec x \vec{\bar x}^T@f$ to a matrix @f$A@f$.
 *
 * The Computation @f$A \gets A + \alpha \vec x \vec{\bar x}^T@f$ is performed.
 *
 * @param uplo Determines if an "upper" right or a "lower" left triangular matrix
 *   @f$A@f$ should be used.
 * @param n Number of rows and columns for the matrix @p a.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param x First field vector @f$\vec x@f$. Dimension of @p x should be @p m.
 * @param incx Stride for elements of @p x.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cher_(const char *uplo, const unsigned *n, const float *alpha,
    const float _Complex *x, const unsigned *incx, float _Complex *a,
    const unsigned *lda);
/** @endcond */
#define h2_syr(uplo, n, alpha, x, incx, a, lda) cher_(uplo, n, alpha, x, incx, a, lda)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
ssyr_(const char *uplo, const unsigned *n, const float *alpha, const float *x,
    const unsigned *incx, float *a, const unsigned *lda);
/** @endcond */
#define h2_syr(uplo, n, alpha, x, incx, a, lda) ssyr_(uplo, n, alpha, x, incx, a, lda)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zher_(const char *uplo, const unsigned *n, const double *alpha,
    const double _Complex *x, const unsigned *incx, double _Complex *a,
    const unsigned *lda);
/** @endcond */
#define h2_syr(uplo, n, alpha, x, incx, a, lda) zher_(uplo, n, alpha, x, incx, a, lda)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dsyr_(const char *uplo, const unsigned *n, const double *alpha, const double *x,
    const unsigned *incx, double *a, const unsigned *lda);
/** @endcond */
#define h2_syr(uplo, n, alpha, x, incx, a, lda) dsyr_(uplo, n, alpha, x, incx, a, lda)
#endif
#endif

/****************************************************
 * BLAS level 3
 ****************************************************/

/****************************************************
 * _GEMM
 ****************************************************/

/**
 * @brief Compute matrix-matrix-multiplication for general matrices @f$A,B@f$
 *   and @f$C@f$.
 *
 *   Depending on the values of @p transa and @p transb one of the following
 *   operations is performed:
 *   @f$C \gets \beta C + \alpha A B@f$ or
 *   @f$C \gets \beta C + \alpha \bar A^T B@f$ or
 *   @f$C \gets \beta C + \alpha A \bar B^T@f$ or
 *   @f$C \gets \beta C + \alpha \bar A^T \bar B^T@f$.
 *
 * @param transa If set to "Conjugate transposed" @f$\bar A^T@f$ will be used.<br>
 *               If set to "Transposed" @f$A^T@f$ will be used.<br>
 *               In all other cases @f$A@f$ will be used for the computation.
 * @param transb If set to "Conjugate transposed" @f$\bar B^T@f$ will be used.<br>
 *               If set to "Transposed" @f$B^T@f$ will be used.<br>
 *               In all other cases @f$B@f$ will be used for the computation.
 * @param m Number of rows for the matrix @p a and matrix @p c.
 * @param n Number of columns for the matrix @p b and the matrix c.
 * @param k Number of columns for the matrix @p a and rows of the matrix @p b.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param b The matrix @f$B@f$.
 * @param ldb Leading dimension of the matrix @f$B@f$
 * @param beta Field scalar value used in the above mentioned computation.
 * @param c The matrix @f$C@f$.
 * @param ldc Leading dimension of the matrix @f$C@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgemm_(const char *transa, const char *transb, const unsigned *m,
    const unsigned *n, const unsigned *k, const float _Complex *alpha,
    const float _Complex *a, const unsigned *lda, const float _Complex *b,
    const unsigned *ldb, const float _Complex *beta, float _Complex *c,
    const unsigned *ldc);
/** @endcond */
#define h2_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) cgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sgemm_(const char *transa, const char *transb, const unsigned *m,
    const unsigned *n, const unsigned *k, const float *alpha, const float *a,
    const unsigned *lda, const float *b, const unsigned *ldb,
    const float *beta, float *c, const unsigned *ldc);
/** @endcond */
#define h2_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgemm_(const char *transa, const char *transb, const unsigned *m,
    const unsigned *n, const unsigned *k, const double _Complex *alpha,
    const double _Complex *a, const unsigned *lda, const double _Complex *b,
    const unsigned *ldb, const double _Complex *beta, double _Complex *c,
    const unsigned *ldc);
/** @endcond */
#define h2_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dgemm_(const char *transa, const char *transb, const unsigned *m,
    const unsigned *n, const unsigned *k, const double *alpha, const double *a,
    const unsigned *lda, const double *b, const unsigned *ldb,
    const double *beta, double *c, const unsigned *ldc);
/** @endcond */
#define h2_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
#endif
#endif

/****************************************************
 * _TRMM
 ****************************************************/

/**
 * @brief Compute matrix-matrix-multiplication for triangular matrices @f$A@f$
 *  and a general matrix @f$B@f$.
 *
 *   Depending on the values of @p transa and @p side one of the following
 *   operations is performed:
 *   @f$B \gets \alpha A B@f$ or
 *   @f$B \gets \alpha \bar A^T B@f$ or
 *   @f$B \gets \alpha B A@f$ or
 *   @f$B \gets \alpha B \bar A^T@f$.
 *
 *   The operation is performed in-place, i.e. matrix @f$B@f$ will be overwritten
 *   by the product of the both matrices.
 *
 * @param side Determines if the matrix @f$A@f$ should be multiplied from "Left"
 *   or from "Right" to the matrix @f$B@f$.
 * @param uplo Determines if an "upper" right or a "lower" left triangular matrix
 *   @f$A@f$ should be used.
 * @param transa If set to "Conjugate transposed" @f$\bar A^T@f$ will be used.<br>
 *               If set to "Transposed" @f$A^T@f$ will be used.<br>
 *               In all other cases @f$A@f$ will be used for the computation.
 * @param diag Determines if a "Unit" diagonal should be used or "Non-unit"
 *        matrix coefficients stored within the matrix @f$A@f$.
 * @param m Number of rows for the matrix @p b.
 * @param n Number of columns for the matrix @p b.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param b The matrix @f$B@f$.
 * @param ldb Leading dimension of the matrix @f$B@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ctrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const float _Complex *alpha,
    const float _Complex *a, const unsigned *lda, float _Complex *b,
    const unsigned *ldb);
/** @endcond */
#define h2_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) ctrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
strmm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const float *alpha, const float *a,
    const unsigned *lda, float *b, const unsigned *ldb);
/** @endcond */
#define h2_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) strmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ztrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const double _Complex *alpha,
    const double _Complex *a, const unsigned *lda, double _Complex *b,
    const unsigned *ldb);
/** @endcond */
#define h2_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) ztrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const double *alpha, const double *a,
    const unsigned *lda, double *b, const unsigned *ldb);
/** @endcond */
#define h2_trmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) dtrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#endif
#endif

/****************************************************
 * _TRSM
 ****************************************************/

/**
 * @brief Solves a linear system of equations with a triangular matrix
 * @f$A@f$ and a right-hand-side matrix @f$B@f$.
 *
 * Depending on the values of @p side and @p transa one of the following systems
 * is being solved for a right-hand-side @f$B@f$ in-place:<br>
 * @f$B \gets \alpha A^{-1} B@f$ or
 * @f$B \gets \alpha \bar A^{-T} B@f$ or
 * @f$B \gets \alpha B A^{-1}@f$ or
 * @f$B \gets \alpha B \bar A^{-T}@f$.
 *
 * @param side Determines if the matrix @f$A^{-1}@f$ should be multiplied from "Left"
 *   or from "Right" to the matrix @f$B@f$.
 * @param uplo Determines if an "upper" right or a "lower" left triangular matrix
 *   @f$A@f$ should be used.
 * @param transa If set to "Conjugate transposed" @f$\bar A^T@f$ will be used.<br>
 *               If set to "Transposed" @f$A^T@f$ will be used.<br>
 *               In all other cases @f$A@f$ will be used for the computation.
 * @param diag Determines if a "Unit" diagonal should be used or "Non-unit"
 *        matrix coefficients stored within the matrix @f$A@f$.
 * @param m Number of rows for the matrix @p b.
 * @param n Number of columns for the matrix @p b.
 * @param alpha Field scalar value used in the above mentioned computation.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param b The matrix @f$B@f$.
 * @param ldb Leading dimension of the matrix @f$B@f$
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ctrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const float _Complex *alpha,
    const float _Complex *a, const unsigned *lda, _Complex float *b,
    const unsigned *ldb);
/** @endcond */
#define h2_trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) ctrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
strsm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const float *alpha, const float *a,
    const unsigned *lda, float *b, const unsigned *ldb);
/** @endcond */
#define h2_trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) strsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
ztrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const double _Complex *alpha,
    const double _Complex *a, const unsigned *lda, _Complex double *b,
    const unsigned *ldb);
/** @endcond */
#define h2_trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) ztrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
    const unsigned *m, const unsigned *n, const double *alpha, const double *a,
    const unsigned *lda, double *b, const unsigned *ldb);
/** @endcond */
#define h2_trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#endif
#endif

/****************************************************
 * LAPACK
 ****************************************************/

/****************************************************
 * _LACGV
 ****************************************************/

/**
 * @brief Conjugates a vector of field @f$\vec x@f$.
 *
 * @f$\vec x \gets \vec \bar{x}@f$.
 *
 * @param n Length of the vector @f$\vec x@f$.
 * @param x Field vector @f$\vec x@f$.
 * @param incx Stride for elements of @p x.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
clacgv_(const unsigned *n, float _Complex *x, const unsigned *incx);
/** @endcond */
#define h2_lacgv(n, x, incx) clacgv_(n, x, incx)
#else
#define h2_lacgv(n, x, incx) (void)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zlacgv_(const unsigned *n, double _Complex *x, const unsigned *incx);
/** @endcond */
#define h2_lacgv(n, x, incx) zlacgv_(n, x, incx)
#else
#define h2_lacgv(n, x, incx) (void)
#endif
#endif

/****************************************************
 * _POTRF
 ****************************************************/

/**
 * @brief Computes the @em Cholesky @em factorization of symmetric/hermitian
 * positive definite matrix @f$A@f$ in-place.
 *
 * The matrix @p a gets partially overwritten by either a lower triangular matrix
 * @f$L@f$ if @p uplo is set to "Lower" or by an upper triangular matrix @f$U@f$,
 * if @p uplo is set to "Upper". In either cases it holds @f$A = L \bar L^T@f$ or
 * @f$A = \bar U^T U@f$.
 *
 * @param uplo Determines whether the lower triangular matrix @f$L@f$ or the
 * upper triangular matrix @f$U@f$ should be stored in the matrix @p a. The
 * value "Lower" refers to the first and the value "Upper" refers to the latter
 * case.
 * @param n Number of rows and columns of the matrix @f$A@f$.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param info Return value for the computation. Equals zero if everything
 * worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cpotrf_(const char *uplo, const unsigned *n, float _Complex *a,
    const unsigned *lda, int *info);
/** @endcond */
#define h2_potrf(uplo, n, a, lda, info) cpotrf_(uplo, n, a, lda, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
spotrf_(const char *uplo, const unsigned *n, float *a, const unsigned *lda,
    int *info);
/** @endcond */
#define h2_potrf(uplo, n, a, lda, info) spotrf_(uplo, n, a, lda, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zpotrf_(const char *uplo, const unsigned *n, double _Complex *a,
    const unsigned *lda, int *info);
/** @endcond */
#define h2_potrf(uplo, n, a, lda, info) zpotrf_(uplo, n, a, lda, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dpotrf_(const char *uplo, const unsigned *n, double *a, const unsigned *lda,
    int *info);
/** @endcond */
#define h2_potrf(uplo, n, a, lda, info) dpotrf_(uplo, n, a, lda, info)
#endif
#endif

/****************************************************
 * _LARF / _LARFG
 ****************************************************/

/**
 * @brief Applies an elementary reflector @f$H = I - \tau \vec v \vec v^T@f$ to
 * a matrix @f$C@f$.
 *
 * Depending on the value of @p side either @f$C \gets H C@f$ or
 * @f$C \gets C H@f$ is being performed.
 *
 * @param side Determines if the matrix @f$H@f$ should be multiplied from "Left"
 *   or from "Right" to the matrix @f$C@f$.
 * @param m Number of rows for the matrix @f$C@f$.
 * @param n Number of columns for the matrix @f$C@f$.
 * @param v Vector of fields @f$\vec v@f$ representing the matrix @f$H@f$.
 * @param incv Stride for the vector @p v.
 * @param tau Field scalar value used in the above mentioned computation of
 *        @f$H@f$.
 * @param c The matrix @f$C@f$.
 * @param ldc Leading dimension of the matrix @f$C@f$
 * @param work Intermediate storage of size @p n if @p side is equal to "Left"
 *    or of size @p m if @p side is equal to "Right".
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
clarf_(const char *side, const unsigned *m, const unsigned *n,
    const float _Complex *v, const unsigned *incv, const float _Complex *tau,
    float _Complex *c, const unsigned *ldc, float _Complex *work);
/** @endcond */
#define h2_larf(side, m, n, v, incv, tau, c, ldc, work) clarf_(side, m, n, v, incv, tau, c, ldc, work)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
slarf_(const char *side, const unsigned *m, const unsigned *n, const float *v,
    const unsigned *incv, const float *tau, float *c, const unsigned *ldc,
    float *work);
/** @endcond */
#define h2_larf(side, m, n, v, incv, tau, c, ldc, work) slarf_(side, m, n, v, incv, tau, c, ldc, work)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zlarf_(const char *side, const unsigned *m, const unsigned *n,
    const double _Complex *v, const unsigned *incv, const double _Complex *tau,
    double _Complex *c, const unsigned *ldc, double _Complex *work);
/** @endcond */
#define h2_larf(side, m, n, v, incv, tau, c, ldc, work) zlarf_(side, m, n, v, incv, tau, c, ldc, work)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dlarf_(const char *side, const unsigned *m, const unsigned *n, const double *v,
    const unsigned *incv, const double *tau, double *c, const unsigned *ldc,
    double *work);
/** @endcond */
#define h2_larf(side, m, n, v, incv, tau, c, ldc, work) dlarf_(side, m, n, v, incv, tau, c, ldc, work)
#endif
#endif

/**
 * @brief Create an elementary reflector @f$H@f$.
 *
 * The properties of @f$H@f$ are defined as
 * @f[\bar H^T \begin{pmatrix} \alpha \\ \vec x \end{pmatrix}
 *  = \begin{pmatrix} \beta \\ \vec 0 \end{pmatrix}, \qquad
 *  \bar H^T H = I. @f]
 *
 * @f$H@f$ will be computed as @f[H = I - \tau
 * \begin{pmatrix} 1 \\ \vec v \end{pmatrix}
 * \begin{pmatrix} 1 & \bar \vec v^T \end{pmatrix}. @f]
 *
 * @param n Dimension of the elementary reflector.
 * @param alpha On entry the scalar field value @f$\alpha@f$, on exit the
 *   scalar field value @f$\beta@f$.
 * @param x On entry the field vector @f$\vec x@f$, on exit the field vector
 *   @f$\vec v@f$.
 * @param incx Stride for the vector @p x.
 * @param tau The field scalar @f$\tau@f$.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
clarfg_(const unsigned *n, float _Complex *alpha, float _Complex *x,
    const unsigned *incx, float _Complex *tau);
/** @endcond */
#define h2_larfg(n, alpha, x, incx, tau) clarfg_(n, alpha, x, incx, tau)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
slarfg_(const unsigned *n, float *alpha, float *x, const unsigned *incx, float *tau);
/** @endcond */
#define h2_larfg(n, alpha, x, incx, tau) slarfg_(n, alpha, x, incx, tau)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zlarfg_(const unsigned *n, double _Complex *alpha, double _Complex *x,
    const unsigned *incx, double _Complex *tau);
/** @endcond */
#define h2_larfg(n, alpha, x, incx, tau) zlarfg_(n, alpha, x, incx, tau)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dlarfg_(const unsigned *n, double *alpha, double *x, const unsigned *incx,
    double *tau);
/** @endcond */
#define h2_larfg(n, alpha, x, incx, tau) dlarfg_(n, alpha, x, incx, tau)
#endif
#endif

/****************************************************
 * _GEQRF
 ****************************************************/

/**
 * @brief Compute a QR-decomposition of the matrix @f$A@f$.
 *
 * @f$A@f$ will be factorized into an orthogonal/unitary matrix @f$Q@f$ and
 * an upper triangular matrix @f$R@f$.<br>
 * Both matrices will be stored into the original matrix @f$A@f$ and an additional
 * vector of field @f$\vec \tau@f$.
 *
 * @param m Rows of the matrix @f$A@f$.
 * @param n Columns of the matrix @f$A@f$.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param tau Field vector containing the scaling factors for the elementary
 *   reflectors.
 * @param work Intermediate storage of size @p lwork.
 * @param lwork Size of the intermediate storage @p work. For a detailed
 *   description of this parameter, please  refer to LAPACK documentation.
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgeqrf_(const unsigned *m, const unsigned *n, float _Complex *a,
    const unsigned *lda, float _Complex *tau, float _Complex *work,
    const int *lwork, int *info);
/** @endcond */
#define h2_geqrf(m, n, a, lda, tau, work, lwork, info) cgeqrf_(m, n, a, lda, tau, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sgeqrf_(const unsigned *m, const unsigned *n, float *a, const unsigned *lda,
    float *tau, float *work, const int *lwork, int *info);
/** @endcond */
#define h2_geqrf(m, n, a, lda, tau, work, lwork, info) sgeqrf_(m, n, a, lda, tau, work, lwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgeqrf_(const unsigned *m, const unsigned *n, double _Complex *a,
    const unsigned *lda, double _Complex *tau, double _Complex *work,
    const int *lwork, int *info);
/** @endcond */
#define h2_geqrf(m, n, a, lda, tau, work, lwork, info) zgeqrf_(m, n, a, lda, tau, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dgeqrf_(const unsigned *m, const unsigned *n, double *a, const unsigned *lda,
    double *tau, double *work, const int *lwork, int *info);
/** @endcond */
#define h2_geqrf(m, n, a, lda, tau, work, lwork, info) dgeqrf_(m, n, a, lda, tau, work, lwork, info)
#endif
#endif

/****************************************************
 * _ORMQR
 ****************************************************/

/**
 * @brief Apply orthogonal/unitary matrix @f$Q@f$ from "Left" or from "Right" to
 * a matrix @f$C@f$.
 *
 * Depending on the values of @p side and @p trans one of the following operations
 * is performed: @f$C \gets Q C@f$ or @f$C \gets \bar Q^T C@f$ or
 * @f$C \gets C Q@f$ or @f$C \gets C \bar Q^T@f$.
 *
 * @f$Q@f$ is represented as product of elementary reflectors
 * @f$H(1) \cdot \ldots \cdot H(k)@f$.
 *
 * @param side Determines if the matrix @f$H@f$ should be multiplied from "Left"
 *   or from "Right" to the matrix @f$C@f$.
 * @param trans If set to "Conjugate transposed" @f$\bar Q^T@f$ will be used.<br>
 *               If set to "Transposed" @f$Q^T@f$ will be used.<br>
 *               In all other cases @f$Q@f$ will be used for the computation.
 * @param m Number of rows for the matrix @p C.
 * @param n Number of columns for the matrix @p C.
 * @param k Number of elementary reflectors.
 * @param a Matrix that contains the elementary reflectors of @f$Q@f$ as returned
 *   by @ref h2_geqrf.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param tau Field vector that contains the scaling factors of the elementary
 *   reflectors.
 * @param c Matrix @f$C@f$ that will be overwritten by the above mentioned
 *   computation.
 * @param ldc Leading dimension of the matrix @f$C@f$
 * @param work Intermediate storage of size @p lwork.
 * @param lwork Size of the intermediate storage @p work. For a detailed
 *   description of this parameter, please  refer to LAPACK documentation.
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cunmqr_(const char *side, const char *trans, const unsigned *m,
    const unsigned *n, const unsigned *k, const float _Complex *a,
    const unsigned *lda, const float _Complex *tau, float _Complex *c,
    const unsigned *ldc, float _Complex *work, const int *lwork, int *info);
/** @endcond */
#define h2_ormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info) cunmqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sormqr_(const char *side, const char *trans, const unsigned *m,
    const unsigned *n, const unsigned *k, const float *a, const unsigned *lda,
    const float *tau, float *c, const unsigned *ldc, float *work,
    const int *lwork, int *info);
/** @endcond */
#define h2_ormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info) sormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zunmqr_(const char *side, const char *trans, const unsigned *m,
    const unsigned *n, const unsigned *k, const double _Complex *a,
    const unsigned *lda, const double _Complex *tau, double _Complex *c,
    const unsigned *ldc, double _Complex *work, const int *lwork, int *info);
/** @endcond */
#define h2_ormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info) zunmqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dormqr_(const char *side, const char *trans, const unsigned *m,
    const unsigned *n, const unsigned *k, const double *a, const unsigned *lda,
    const double *tau, double *c, const unsigned *ldc, double *work,
    const int *lwork, int *info);
/** @endcond */
#define h2_ormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info) dormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
#endif
#endif

/****************************************************
 * _ORGQR
 ****************************************************/

/**
 * @brief Generate the orthogonal/unitary matrix @f$Q@f$ out of its elementary
 *   reflectors @f$Q = H(1) \cdot \ldots \cdot H(k)@f$.
 *
 * The matrix @f$Q@f$ will overwrite the matrix @f$A@f$ on exit.
 *
 * @param m Number of rows for the matrix @p C.
 * @param n Number of columns for the matrix @p C.
 * @param k Number of elementary reflectors.
 * @param a Matrix that contains the elementary reflectors of @f$Q@f$ as returned
 *   by @ref h2_geqrf on entry and will contain the matrix @f$Q@f$ at exit.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param tau Field vector that contains the scaling factors of the elementary
 *   reflectors.
 * @param work Intermediate storage of size @p lwork.
 * @param lwork Size of the intermediate storage @p work. For a detailed
 *   description of this parameter, please  refer to LAPACK documentation.
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cungqr_(const unsigned *m, const unsigned *n, const unsigned *k,
    float _Complex *a, const unsigned *lda, float _Complex *tau,
    float _Complex *work, const unsigned *lwork, unsigned *info);
/** @endcond */
#define h2_orgqr(m, n, k, a, lda, tau, work, lwork, info) cungqr_(m, n, k, a, lda, tau, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sorgqr_(const unsigned *m, const unsigned *n, const unsigned *k, float *a,
    const unsigned *lda, float *tau, float *work, const unsigned *lwork,
    unsigned *info);
/** @endcond */
#define h2_orgqr(m, n, k, a, lda, tau, work, lwork, info) sorgqr_(m, n, k, a, lda, tau, work, lwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zungqr_(const unsigned *m, const unsigned *n, const unsigned *k,
    double _Complex *a, const unsigned *lda, double _Complex *tau,
    double _Complex *work, const unsigned *lwork, unsigned *info);
/** @endcond */
#define h2_orgqr(m, n, k, a, lda, tau, work, lwork, info) zungqr_(m, n, k, a, lda, tau, work, lwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dorgqr_(const unsigned *m, const unsigned *n, const unsigned *k, double *a,
    const unsigned *lda, double *tau, double *work, const unsigned *lwork,
    unsigned *info);
/** @endcond */
#define h2_orgqr(m, n, k, a, lda, tau, work, lwork, info) dorgqr_(m, n, k, a, lda, tau, work, lwork, info)
#endif
#endif

/****************************************************
 * _STEQR
 ****************************************************/

/**
 * @brief Compute all eigenvalues and optionally all eigenvectors of a
 * symmetric tridiagonal matrix @f$A@f$ using implicit QL or QR method.
 *
 * @param compz Is equal to "No Vectors" if only the eigenvalues should be computed
 *   and equal to "Vectors" if also eigenvectors should be computed. Otherwise
 *   @p compz should be equal to "Identity" if @f$Z@f$ was initialized by an
 *   identity matrix and both eigenvalues and eigenvectors should be computed.
 * @param n Order of the matrix @f$A@f$.
 * @param d Coefficients for the diagonal of @f$A@f$.
 * @param e Coefficients for the subdiagonals of @f$A@f$.
 * @param z If Eigenvectors should be computed, the columns of the matrix
 *   @f$Z@f$ will contain them at exit.
 * @param ldz Leading dimension of the matrix @f$Z@f$.
 * @param work Intermediate storage of size @p lwork.
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
csteqr_(const char *compz, const unsigned *n, float *d, float *e,
    float _Complex *z, const unsigned *ldz, float *work, int *info);
/** @endcond */
#define h2_steqr(compz, n, d, e, z, ldz, work, info) csteqr_(compz, n, d, e, z, ldz, work, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
ssteqr_(const char *compz,const unsigned *n,float *d,float *e, float *z,
    const unsigned *ldz, float *work, int *info);
/** @endcond */
#define h2_steqr(compz, n, d, e, z, ldz, work, info) ssteqr_(compz, n, d, e, z, ldz, work, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zsteqr_(const char *compz, const unsigned *n, double *d, double *e,
    double _Complex *z, const unsigned *ldz, double *work, int *info);
/** @endcond */
#define h2_steqr(compz, n, d, e, z, ldz, work, info) zsteqr_(compz, n, d, e, z, ldz, work, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dsteqr_(const char *compz, const unsigned *n, double *d, double *e, double *z,
    const unsigned *ldz, double *work, int *info);
/** @endcond */
#define h2_steqr(compz, n, d, e, z, ldz, work, info) dsteqr_(compz, n, d, e, z, ldz, work, info)
#endif
#endif

/****************************************************
 * _STEV
 ****************************************************/

/**
 * @brief Compute all eigenvalues and optionally all eigenvectors of a real
 * symmetric tridiagonal matrix @f$A@f$.
 *
 * @param jobz Is equal to "No Vectors" if only the eigenvalues should be computed
 *   and equal to "Vectors" if also eigenvectors should be computed.
 * @param n Order of the matrix @f$A@f$.
 * @param d Coefficients for the diagonal of @f$A@f$.
 * @param e Coefficients for the subdiagonals of @f$A@f$.
 * @param z If Eigenvectors should be computed, the columns of the matrix
 *   @f$Z@f$ will contain them at exit.
 * @param ldz Leading dimension of the matrix @f$Z@f$.
 * @param work Intermediate storage of size @p lwork.
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
/** @cond IMPORT */
IMPORT_PREFIX void
sstev_(const char *jobz, const unsigned *n, float *d, float *e, float *z,
    const unsigned *ldz, float *work, int *info);
/** @endcond */
#define h2_stev(jobz, n, d, e, z, ldz, work, info) sstev_(jobz, n, d, e, z, ldz, work, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dstev_(const char *jobz, const unsigned *n, double *d, double *e, double *z,
    const unsigned *ldz, double *work, int *info);
/** @endcond */
#define h2_stev(jobz, n, d, e, z, ldz, work, info) dstev_(jobz, n, d, e, z, ldz, work, info)
#endif

/****************************************************
 * _GESVD
 ****************************************************/

/**
 * @brief Compute the SVD of a general matrix @f$A@f$
 *
 * @f$A@f$ is decomposed as @f$A = U \Sigma \bar V^T@f$.
 *
 * @param jobu If set to "All" all @p m columns of @f$U@f$ are computed.<br>
 *   If set to "Skinny" only the first min(@p m, @p n) columns of @f$U@f$ are computed.<br>
 *   If set to "Overwrite" only the first min(@p m, @p n) columns of @f$U@f$ are
 *   computed and overwrite the entries of @f$A@f$.<br>
 *   If set to "No Vectors" no left singular vectors will be computed.
 * @param jobvt If set to "All" all @p n columns of @f$\bar V^T@f$ are computed.<br>
 *   If set to "Skinny" only the first min(@p m, @p n) columns of @f$\bar V^T@f$ are computed.<br>
 *   If set to "Overwrite" only the first min(@p m, @p n) columns of @f$\bar V^T@f$ are
 *   computed and overwrite the entries of @f$A@f$.<br>
 *   If set to "No Vectors" no right singular vectors will be computed.
 * @param m Number of rows of the matrix @f$A@f$.
 * @param n Number of rows of the matrix @f$A@f$.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of @f$A@f$.
 * @param s Vector of reals containing the singular vectors on exit.
 * @param u Matrix containing the left singular vectors depending on the choice
 *   of @p jobu on exit.
 * @param ldu Leading dimension of @f$U@f$.
 * @param vt Matrix containing the right singular vectors depending on the choice
 *   of @p jobvt on exit.
 * @param ldvt Leading dimension of @f$\bar V^T@f$.
 * @param work Intermediate storage of size @p lwork.
 * @param lwork Size of the intermediate storage @p work. For a detailed
 *   description of this parameter, please  refer to LAPACK documentation.
 * @param rwork Intermediate storage of size 5*min(@p m, @p n).
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cgesvd_(const char *jobu, const char *jobvt, const unsigned *m,
    const unsigned *n, float _Complex *a, const unsigned *lda, float *s,
    float _Complex *u, const unsigned *ldu, float _Complex *vt,
    const unsigned *ldvt, float _Complex *work, const unsigned *lwork,
    float *rwork, int *info);
/** @endcond */
#define h2_gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info) cgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sgesvd_(const char *jobu, const char *jobvt, const unsigned *m,
    const unsigned *n, float *a, const unsigned *lda, float *s, float *u,
    const unsigned *ldu, float *vt, const unsigned *ldvt, float *work,
    const unsigned *lwork, int *info);
/** @endcond */
#define h2_gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info) sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zgesvd_(const char *jobu, const char *jobvt, const unsigned *m,
    const unsigned *n, double _Complex *a, const unsigned *lda, double *s,
    double _Complex *u, const unsigned *ldu, double _Complex *vt,
    const unsigned *ldvt, double _Complex *work, const unsigned *lwork,
    double *rwork, int *info);
/** @endcond */
#define h2_gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info) zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dgesvd_(const char *jobu, const char *jobvt, const unsigned *m,
    const unsigned *n, double *a, const unsigned *lda, double *s, double *u,
    const unsigned *ldu, double *vt, const unsigned *ldvt, double *work,
    const unsigned *lwork, int *info);
/** @endcond */
#define h2_gesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info) dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
#endif
#endif

/****************************************************
 * _DBSQR
 ****************************************************/

/**
 * @brief Compute the SVD of a real "Lower" or "Upper" bidiagonal matrix @f$B@f$.
 *
 * @f$B@f$ is decomposed as @f$B = Q \Sigma \bar P^T@f$.
 *
 * If unitary matrices @f$U@f$ ad @f$\bar V^T@f$ that reduce a matrix @f$A@f$ to
 * bidiagonal form, then the SVD of @f$A@f$ is given by
 * @f$A = (U Q) \Sigma (\bar P^T \bar V^T)@f$.
 *
 * Optionally the function may also compute @f$\bar Q^T C@f$ for some input matrix
 * @f$C@f$.
 *
 * @param uplo Determines whether an "Upper" or a "Lower" bidiagonal matrix
 *   @f$B@f$ is given.
 * @param n The order of the matrix @f$B@f$.
 * @param ncvt Number of columns of the matrix @f$\bar V^T@f$.
 * @param nru Number of row of the matrix @f$U@f$.
 * @param ncc Number of columns of the matrix @f$C@f$.
 * @param d Real vector of length @p n containing the diagonal entries of @f$B@f$
 *   on entry and the singular values in descending order on exit.
 * @param e Real vector of length @p n - 1 containing the offdiagonal entries of
 *   @f$B@f$.
 * @param vt Matrix of fields @f$\bar V^T@f$ on entry, on exit the matrix
 *   @f$\bar P^T \bar V^T@f$.
 * @param ldvt Leading dimension of the matrix @f$\bar V^T@f$.
 * @param u Matrix of fields @f$U@f$ on entry, on exit the matrix
 *   @f$U Q@f$.
 * @param ldu Leading dimension of the matrix @f$U@f$.
 * @param c Matrix of fields @f$C@f$ on entry, on exit the matrix
 *   @f$\bar Q^T C@f$.
 * @param ldc Leading dimension of the matrix @f$C@f$.
 * @param rwork Intermediate storage of size 2 @p n or max(1, 4*@p n - 4).
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cbdsqr_(const char *uplo, const unsigned *n, const unsigned *ncvt,
    const unsigned *nru, const unsigned *ncc, float *d, float *e,
    float _Complex *vt, const unsigned *ldvt, float _Complex *u,
    const unsigned *ldu, float _Complex *c, const unsigned *ldc, float *rwork,
    int *info);
/** @endcond */
#define h2_bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info) cbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
sbdsqr_(const char *uplo, const unsigned *n, const unsigned *ncvt,
    const unsigned *nru, const unsigned *ncc, float *d, float *e, float *vt,
    const unsigned *ldvt, float *u, const unsigned *ldu, float *c,
    const unsigned *ldc, float *rwork, int *info);
/** @endcond */
#define h2_bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info) sbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zbdsqr_(const char *uplo, const unsigned *n, const unsigned *ncvt,
    const unsigned *nru, const unsigned *ncc, double *d, double *e,
    double _Complex *vt, const unsigned *ldvt, double _Complex *u,
    const unsigned *ldu, double _Complex *c, const unsigned *ldc, double *rwork,
    int *info);
/** @endcond */
#define h2_bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info) zbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dbdsqr_(const char *uplo, const unsigned *n, const unsigned *ncvt,
    const unsigned *nru, const unsigned *ncc, double *d, double *e, double *vt,
    const unsigned *ldvt, double *u, const unsigned *ldu, double *c,
    const unsigned *ldc, double *rwork, int *info);
/** @endcond */
#define h2_bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info) dbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
#endif
#endif

/****************************************************
 * _SYEV / _HEEV
 ****************************************************/

/**
 * @brief Compute all eigenvalues and optionally all eigenvectors of a
 * hermitian / symmetric matrix @f$A@f$.
 *
 * @param jobz Is equal to "No Vectors" if only the eigenvalues should be computed
 *   and equal to "Vectors" if also eigenvectors should be computed.
 * @param uplo Determines whether the lower triangular or the
 * upper triangular part of the matrix @f$A@f$ should be used.
 * The value "Lower" refers to the first and the value "Upper" refers to
 * the latter case.
 * @param n Number of rows and columns of the matrix @f$A@f$.
 * @param a The matrix @f$A@f$.
 * @param lda Leading dimension of the matrix @f$A@f$
 * @param w On exit the field vector @f$w@f$ contains the eigenvalues in
 *   ascending order.
 * @param work Intermediate storage of size @p lwork.
 * @param lwork Size of the intermediate storage @p work. For a detailed
 *   description of this parameter, please  refer to LAPACK documentation.
 * @param rwork Intermediate storage of size max(1, @p 3*n-2).
 * @param info Return value for the computation. Equals zero if everything
 *  worked fine, otherwise a return code that indicates the problem is returned.
 */
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
cheev_(const char *jobz, const char *uplo, const unsigned *n, float _Complex *a,
    const unsigned *lda, const float *w, float _Complex *work,
    const unsigned *lwork, float *rwork, int *info);
/** @endcond */
#define h2_heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) cheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
ssyev_(const char *jobz, const char *uplo, const unsigned *n, float *a,
    const unsigned *lda, const float *w, float *work, const unsigned *lwork,
    int *info);
/** @endcond */
#define h2_heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info)
#endif
#else
#ifdef USE_COMPLEX
/** @cond IMPORT */
IMPORT_PREFIX void
zheev_(const char *jobz, const char *uplo, const unsigned *n,
    double _Complex *a, const unsigned *lda, const double *w,
    double _Complex *work, const unsigned *lwork, double *rwork, int *info);
/** @endcond */
#define h2_heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
#else
/** @cond IMPORT */
IMPORT_PREFIX void
dsyev_(const char *jobz, const char *uplo, const unsigned *n, double *a,
    const unsigned *lda, const double *w, double *work, const unsigned *lwork,
    int *info);
/** @endcond */
#define h2_heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info) dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info)
#endif
#endif

#endif

/**
 * @}
 */

#endif

