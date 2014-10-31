/* ------------------------------------------------------------
 This is the file "eigensolvers.h" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

/**
 * @file eigensolvers.h
 * @author Steffen B&ouml;rm
 */

#ifndef EIGENSOLVERS_H
#define EIGENSOLVERS_H

/**
 *  @defgroup eigensolvers eigensolvers
 *  @brief Solve symmetric eigenproblems and compute singular value
 *     decompositions.
 *  @{
 */

/** @brief Tridiagonal matrix. */
typedef struct _tridiag tridiag;

/** @brief Pointer to a @ref tridiag object. */
typedef tridiag *ptridiag;

/** @brief Pointer to a constant @ref tridiag object. */
typedef const tridiag *pctridiag;

#include "amatrix.h"
#include "settings.h"

/**
 * @brief Tridiagonal matrix @f$T@f$, represented by vectors containing
 *    the diagonal, sub- and superdiagonal.
 */
struct _tridiag {
  /**
   * @brief Diagonal represented by the vector
   *    @f$(t_{11},t_{22},\ldots,t_{nn})@f$,
   *    dimension <tt>dim</tt>.
   */
  pfield d;

  /**
   * @brief Subdiagonal represented by the vector
   *    @f$(t_{21},t_{32},\ldots,t_{n,n-1})@f$,
   *    dimension <tt>dim-1</tt>.
   */
  pfield l;

  /**
   * @brief Superdiagonal represented by the vector
   *    @f$(t_{12},t_{23},\ldots,t_{n-1,n})@f$,
   *    dimension <tt>dim-1</tt>.
   */
  pfield u;

  /**
   * @brief Matrix dimension.
   */
  uint dim;

  /**
   * @brief If this a submatrix, this points to the supermatrix it
   *    was taken from.
   */
  ptridiag owner;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief Initialize a @ref tridiag object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient vectors.
 *
 *  @remark Should always be matched by a call to @ref uninit_tridiag.
 *
 *  @param T Object to be initialized.
 *  @param dim Dimension of the new matrix.
 *  @returns Initialized @ref tridiag object. */
HEADER_PREFIX ptridiag
init_tridiag(ptridiag T, uint dim);

/**
 * @brief Initialize a @ref tridiag object to represent a submatrix.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref amatrix for the coefficient vectors, leading to a
 *  new matrix representing a submatrix of the source.
 *
 *  @remark Should always be matched by a call to @ref uninit_tridiag
 *  that will <em>not</em> release the coefficient storage.
 *
 *  @param T Object to be initialized.
 *  @param src Source matrix @f$S@f$.
 *  @param dim Dimension of submatrix.
 *  @param off Row and column offset.
 *  @returns Initialized @ref tridiag object. */
HEADER_PREFIX ptridiag
init_sub_tridiag(ptridiag T, ptridiag src, uint dim, uint off);

/** @brief Initialize a @ref tridiag object using pre-allocated
 *  storage.
 *
 *  Sets up the components of the object and uses the storage of
 *  a vector of at least dimension <tt>3*dim-2</tt> for the
 *  coefficient vectors.
 *
 *  @remark Should always be matched by a call to @ref uninit_tridiag
 *  that will <em>not</em> release the coefficient storage.
 *
 *  @param T Object to be initialized.
 *  @param src Source vector.
 *  @param dim Dimension of the new matrix.
 *  @returns Initialized @ref tridiag object. */
HEADER_PREFIX ptridiag
init_vec_tridiag(ptridiag T, pavector src, uint dim);

/** @brief Uninitialize a @ref tridiag object.
 *
 *  Invalidates pointers, freeing corresponding storage if <tt>owner==0</tt>,
 *  and prepares the object for deletion.
 *
 *  @param T Object to be uninitialized. */
HEADER_PREFIX void
uninit_tridiag(ptridiag T);

/** @brief Create a new @ref tridiag object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_tridiag.
 *
 *  @param dim Dimension of the new matrix.
 *  @returns New @ref tridiag object. */
HEADER_PREFIX ptridiag
new_tridiag(uint dim);

/** @brief Delete a @ref tridiag object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they will otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param T Object to be deleted. */
HEADER_PREFIX void
del_tridiag(ptridiag T);

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

/** @brief Copy a matrix into another matrix, @f$T_{\rm copy} \gets T@f$.
 *
 *  @param T Source matrix @f$T@f$.
 *  @param Tcopy Target matrix @f$T_{\rm copy}@f$. */
HEADER_PREFIX void
copy_tridiag(pctridiag T, ptridiag Tcopy);

/** @brief Compute the Frobenius norm @f$\|T - T_s\|_F@f$ of the difference
 *  of a tridiagonal matrix @f$T@f$ and a general matrix @f$T_s@f$.
 *
 *  @param T First matrix, represented as @ref tridiag object.
 *  @param Ts Second matrix, represented as @ref amatrix object.
 *  @returns Frobenius norm @f$\|T - T_s\|_F@f$. */
HEADER_PREFIX real
check_tridiag(pctridiag T, pcamatrix Ts);

/** @brief Compute the Frobenius norm @f$\|T - T_s\|_F@f$ of the difference
 *  of a lower bidiagonal matrix @f$T@f$ and a general matrix @f$T_s@f$.
 *
 *  @param T First matrix, represented as @ref tridiag object, where
 *    the superdiagonal <tt>T->u</tt> is ignored.
 *  @param Ts Second matrix, represented as @ref amatrix object.
 *  @returns Frobenius norm @f$\|T - T_s\|_F@f$. */
HEADER_PREFIX real
check_lower_tridiag(pctridiag T, pcamatrix Ts);

/* ------------------------------------------------------------
 Compute eigenvalues and eigenvectors,
 only for self-adjoint tridiagonal matrices
 ------------------------------------------------------------ */

/** @brief Perform one implicit QR step on a self-adjoint tridiagonal
 *  matrix.
 *
 *  @param T Matrix @f$T@f$, will be overwritten by
 *    @f$\widehat{Q}^* T \widehat{Q}@f$, where @f$\widehat{Q}@f$ is the
 *    unitary transformation corresponding to the QR step.
 *  @param shift Shift parameter @f$\mu@f$, the matrix @f$T@f$ is
 *    shifted to @f$T-\mu I@f$ before computing the QR factorization.
 *  @param Q Can be used to accumulate transformations, @f$Q@f$ will be
 *    overwritten by @f$Q \widehat{Q}@f$.
 *    If <tt>Q==0</tt>, the transformations will not be accumulated. */
HEADER_PREFIX void
qrstep_tridiag(ptridiag T, field shift, pamatrix Q);

/** @brief Solve a self-adjoint tridiagonal eigenproblem.
 *
 *  Eigenvalues will be on the diagonal of <tt>T</tt>, in ascending order.
 *
 *  @remark This is a very simple implementation using Wilkinson shifts.
 *    It is only intended as a stand-in on systems that do not offer
 *    a LAPACK library.
 *    Whenever possible, @ref muleig_tridiag should be used instead.
 * 
 *  @param T Matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = \widehat{Q} D \widehat{Q}^*@f$ with
 *    a unitary matrix @f$\widehat{Q}@f$.
 *  @param Q Unitary matrix @f$Q@f$ used to accumulate transformations,
 *    will be overwritten by @f$Q \widehat{Q}@f$.
 *    If <tt>Q==0</tt>, the transformations will not be accumulated.
 *  @param maxiter Upper bound for the number of QR steps.
 *  @returns Number of QR steps. */
HEADER_PREFIX uint
sb_muleig_tridiag(ptridiag T, pamatrix Q, uint maxiter);

/** @brief Solve a self-adjoint tridiagonal eigenproblem.
 *
 *  Eigenvalues will be on the diagonal of <tt>T</tt>, in ascending order.
 *
 *  @param T Matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = \widehat{Q} D \widehat{Q}^*@f$ with
 *    a unitary matrix @f$\widehat{Q}@f$.
 *  @param Q Unitary matrix @f$Q@f$ used to accumulate transformations,
 *    will be overwritten by @f$Q \widehat{Q}@f$.
 *    If <tt>Q==0</tt>, the transformations will not be accumulated.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
muleig_tridiag(ptridiag T, pamatrix Q);

/** @brief Solve a self-adjoint tridiagonal eigenproblem.
 *
 *  Eigenvalues will be on the diagonal of <tt>T</tt>, in ascending order.
 *
 *  @param T Matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = Q D Q^*@f$ with
 *    a unitary matrix @f$Q@f$.
 *  @param Q If <tt>Q!=0</tt>, this matrix will be filled with an orthonormal
 *    basis of eigenvectors.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
eig_tridiag(ptridiag T, pamatrix Q);

/** @brief Hessenberg tridiagonalization of a self-adjoint matrix.
 *
 *  Compute a unitary matrix @f$Q@f$ and a tridiagonal self-adjoint matrix
 *  @f$T@f$ such that @f$A = Q T Q^*@f$.
 *
 *  @remark This function is only intended as an intermediate step of
 *    @ref eig_amatrix on systems that do not offer a LAPACK library.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param T Tridiagonal matrix @f$T@f$.
 *  @param Q If <tt>Q!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$Q@f$. */
HEADER_PREFIX void
sb_tridiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix Q);

/** @brief Solve a self-adjoint eigenproblem.
 *
 *  @remark This is a very simple implementation using Wilkinson shifts.
 *    It is only intended as a stand-in on systems that do not offer
 *    a LAPACK library.
 *    Whenever possible, @ref eig_amatrix should be used instead.
 * 
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param lambda Eigenvalues of @f$A@f$, in ascending order.
 *  @param Q If <tt>Q!=0</tt>, the columns of this matrix will contain
 *    an orthonormal basis of eigenvectors corresponding to the
 *    eigenvalues in <tt>lambda</tt>.
 *  @param maxiter Upper bound for the number of QR steps.
 *  @returns Number of QR steps. */
HEADER_PREFIX uint
sb_eig_amatrix(pamatrix A, pavector lambda, pamatrix Q, uint maxiter);

/** @brief Solve a self-adjoint eigenproblem, @f$A e = \lambda e@f$.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param lambda Eigenvalues of @f$A@f$, in ascending order.
 *  @param Q If <tt>Q!=0</tt>, the columns of this matrix will contain
 *    an orthonormal basis of eigenvectors corresponding to the
 *    eigenvalues in <tt>lambda</tt>.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
eig_amatrix(pamatrix A, pavector lambda, pamatrix Q);

/** @brief Solve self-adjoint generalized eigenproblem,
 *  @f$A e = \lambda M e@f$.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param M Matrix @f$M@f$, will be overwritten by the Cholesky
 *    factorization of @f$M@f$.
 *  @param lambda Eigenvalues of @f$A@f$, in ascending order.
 *  @param Q If <tt>Q!=0</tt>, the columns of this matrix will contain
 *    an @f$M@f$-orthonormal basis of eigenvectors corresponding to the
 *    eigenvalues in <tt>lambda</tt>.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
geig_amatrix(pamatrix A, pamatrix M, pavector lambda, pamatrix Q);

/* ------------------------------------------------------------
 Singular value decomposition of a sub-bidiagonal matrix
 ------------------------------------------------------------ */

/** @brief Perform one step of the Golub-Kahan iteration on a lower
 *  bidiagonal matrix.
 *
 *  @param T Bidiagonal matrix @f$T@f$, will be overwritten
 *    by @f$\widehat{U}^* T \widehat{V}@f$.
 *  @param shift Shift parameter @f$\mu@f$, the matrix @f$T T^*@f$
 *    is shifted to @f$T T^* - \mu I@f$ before computing the QR factorization.
 *  @param U Can be used to accumulate transformations,
 *    @f$U@f$ will be overwritten by @f$U \widehat{U}@f$.
 *    If <tt>U==0</tt>, the transformations will not be accumulated.
 *  @param Vt Can be used to accumulate transformations,
 *    @f$V_t@f$ will be overwritten by @f$\widehat{V}^* V_t@f$.
 *    If <tt>Vt==0</tt>, the transformations will not be accumulated. */
HEADER_PREFIX void
svdstep_tridiag(ptridiag T, field shift, pamatrix U, pamatrix Vt);

/** @brief Compute the SVD of a bidiagonal matrix.
 *
 *  Singular values will be on the diagonal of <tt>T</tt>, in descending order.
 *
 *  @remark This is a very simple implementation using Wilkinson shifts.
 *    It is only intended as a stand-in on systems that do not offer
 *    a LAPACK library.
 *    Whenever possible, @ref mulsvd_tridiag should be used instead.
 * 
 *  @param T Bidiagonal matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = \widehat{U} D \widehat{V}^*@f$ with
 *    unitary matrices @f$\widehat{U}@f$ and @f$\widehat{V}@f$.
 *  @param U Unitary matrix @f$U@f$ used to accumulate transformations,
 *    will be overwritten by @f$U \widehat{U}@f$.
 *    If <tt>U==0</tt>, the transformations will not be accumulated.
 *  @param Vt Adjoint @f$V_t=V^*@f$ of unitary matrix @f$V@f$ used to
 *    accumulate transformations,
 *    will be overwritten by @f$\widehat{V}^* V_t = (V \widehat{V})^*@f$.
 *    If <tt>Vt==0</tt>, the transformations will not be accumulated.
 *  @param maxiter Upper bound for the number of Golub-Kahan steps.
 *  @returns Number of Golub-Kahan steps. */
HEADER_PREFIX uint
sb_mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt, uint maxiter);

/** @brief Compute the SVD of a bidiagonal matrix.
 *
 *  Compute the factorization @f$T = U \Sigma V^*@f$ with a real diagonal
 *  matrix @f$\Sigma@f$ and unitary matrices @f$U@f$ and @f$V@f$.
 *
 *  @param T Bidiagonal matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = \widehat{U} D \widehat{V}^*@f$ with
 *    unitary matrices @f$\widehat{U}@f$ and @f$\widehat{V}@f$.
 *  @param U Unitary matrix @f$U@f$ used to accumulate transformations,
 *    will be overwritten by @f$U \widehat{U}@f$.
 *    If <tt>U==0</tt>, the transformations will not be accumulated.
 *  @param Vt Adjoint @f$V_t=V^*@f$ of unitary matrix @f$V@f$ used to
 *    accumulate transformations,
 *    will be overwritten by @f$\widehat{V}^* V_t = (V \widehat{V})^*@f$.
 *    If <tt>Vt==0</tt>, the transformations will not be accumulated.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt);

/** @brief Compute the SVD of a bidiagonal matrix.
 *
 *  Singular values will be on the diagonal of <tt>T</tt>, in descending order.
 *
 *  @param T Bidiagonal matrix @f$T@f$, will be overwritten by a diagonal
 *    matrix @f$D@f$ such that @f$T = \widehat{U} D \widehat{V}^*@f$ with
 *    unitary matrices @f$\widehat{U}@f$ and @f$\widehat{V}@f$.
 *  @param U If <tt>U!=0</tt>, the columns will filled by the leading
 *    left singular vectors.
 *  @param Vt If <tt>Vt!=0</tt>, the rows will be filled by the conjugated
 *    leading right singular vectors.
 *  @returns Zero on successful completion, non-zero if the QR iteration
 *    did not converge. */
HEADER_PREFIX uint
svd_tridiag(ptridiag T, pamatrix U, pamatrix Vt);

/** @brief Golub-Kahan bidiagonalization of a matrix.
 *
 *  Compute unitary matrices @f$U@f$ and @f$V@f$ and a lower bidiagonal
 *  matrix @f$T@f$ such that @f$A = U T V^*@f$.
 *
 *  @remark This function is only intended as an intermediate step of
 *    @ref svd_amatrix on systems that do not offer a LAPACK library.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param T Bidiagonal matrix @f$T@f$.
 *  @param U If <tt>U!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$U@f$.
 *  @param Vt If <tt>Vt!=0</tt>, this matrix will be filled with the
 *    adjoint unitary transformation @f$V^*@f$. */
HEADER_PREFIX void
sb_bidiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix U, pamatrix Vt);

/** @brief Compute the SVD of a matrix.
 *
 *  Compute unitary matrices @f$U@f$ and @f$V@f$ and a lower bidiagonal
 *  matrix @f$T@f$ such that @f$A = U T V^*@f$.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param work Auxiliary storage, dimension at least as large as
 *    the maximum of the number of rows and of columns.
 *  @param T Bidiagonal matrix @f$T@f$.
 *  @param U If <tt>U!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$U@f$.
 *  @param Vt If <tt>Vt!=0</tt>, this matrix will be filled with the
 *    adjoint unitary transformation @f$V^*@f$. */
HEADER_PREFIX void
bidiagonalize_amatrix(pamatrix A, pavector work, ptridiag T, pamatrix U,
    pamatrix Vt);

/** @brief Compute the SVD of a matrix.
 *
 *  Compute the factorization @f$A = U \Sigma V^*@f$ with a real diagonal
 *  matrix @f$\Sigma=\mathop{\operatorname{diag}}(\sigma_1,\ldots,\sigma_k)@f$,
 *  @f$\sigma_1\geq\sigma_2\geq\ldots@f$ and unitary matrices @f$U@f$ and
 *  @f$V@f$.
 *
 *  @remark This function is only intended as an intermediate step of
 *    @ref svd_amatrix on systems that do not offer a LAPACK library.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param sigma Singular values @f$\sigma_1,\sigma_2,\ldots@f$.
 *  @param U If <tt>U!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$U@f$.
 *  @param Vt If <tt>Vt!=0</tt>, this matrix will be filled with the
 *    adjoint unitary transformation @f$V^*@f$.
 *  @param maxiter Upper bound for the number of Golub-Kahan steps.
 *  @returns Number of Golub-Kahan steps. */
HEADER_PREFIX uint
sb_svd_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt,
    uint maxiter);

/** @brief Compute the SVD of a matrix.
 *
 *  Compute the factorization @f$A = U \Sigma V^*@f$ with a real diagonal
 *  matrix @f$\Sigma=\mathop{\operatorname{diag}}(\sigma_1,\ldots,\sigma_k)@f$,
 *  @f$\sigma_1\geq\sigma_2\geq\ldots@f$ and unitary matrices @f$U@f$ and
 *  @f$V@f$.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param sigma Singular values @f$\sigma_1,\sigma_2,\ldots@f$.
 *  @param U If <tt>U!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$U@f$.
 *  @param Vt If <tt>Vt!=0</tt>, this matrix will be filled with the
 *    adjoint unitary transformation @f$V^*@f$.
 *  @returns Zero on successful completion, non-zero if the Golub-Kahan
 *    iteration did not converge. */
HEADER_PREFIX uint
svd_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt);

/** @brief Compute the SVD of a matrix and explicitly check the result.
 *
 *  Compute the factorization @f$A = U \Sigma V^*@f$ with a real diagonal
 *  matrix @f$\Sigma=\mathop{\operatorname{diag}}(\sigma_1,\ldots,\sigma_k)@f$,
 *  @f$\sigma_1\geq\sigma_2\geq\ldots@f$ and unitary matrices @f$U@f$ and
 *  @f$V@f$.
 *
 *  @param A Matrix @f$A@f$, will be overwritten by the function.
 *  @param sigma Singular values @f$\sigma_1,\sigma_2,\ldots@f$.
 *  @param U If <tt>U!=0</tt>, this matrix will be filled with the
 *    unitary transformation @f$U@f$.
 *  @param Vt If <tt>Vt!=0</tt>, this matrix will be filled with the
 *    adjoint unitary transformation @f$V^*@f$.
 *  @returns Zero on successful completion, non-zero if the Golub-Kahan
 *    iteration did not converge. */
HEADER_PREFIX uint
svd_verified_amatrix(pamatrix A, pavector sigma, pamatrix U, pamatrix Vt);

/** @} */

#endif
