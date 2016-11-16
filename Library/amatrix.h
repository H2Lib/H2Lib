
/* ------------------------------------------------------------
 * This is the file "amatrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file amatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef AMATRIX_H
#define AMATRIX_H

/** @defgroup amatrix amatrix
 *  @brief Representation of a matrix as an array in column-major order.
 *
 *  The @ref amatrix class is used to handle standard linear algebra
 *  operations like matrix multiplication, factorization, solving
 *  linear systems or eigenvalue problems.
 *  @{ */

/** @brief Representation of a matrix as a column-order array. */
typedef struct _amatrix amatrix;

/** @brief Pointer to @ref amatrix object. */
typedef amatrix *pamatrix;

/** @brief Pointer to constant @ref amatrix object. */
typedef const amatrix *pcamatrix;

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "basic.h"
#include "settings.h"
#include "blas.h"
#include "avector.h"
#include "realavector.h"
#include "krylov.h"

/** @brief Representation of a matrix as an array in column-major order. */
struct _amatrix {
  /** @brief Matrix coefficients in column-major order, i.e., @f$a_{ij}@f$ corresponds to `a[i+j*ld]`.  */
  field *a;

  /** @brief Leading dimension, i.e., increment used to switch from one column to the next.  */
  uint ld;

  /** @brief Number of rows. */
  uint rows;
  /** @brief Number of columns.  */
  uint cols;

  /** @brief Points to owner of coefficient storage if this is a submatrix. */
  void *owner;
};

/* ------------------------------------------------------------ *
 * Constructors and destructors                                 *
 * ------------------------------------------------------------ */

/** @brief Initialize an @ref amatrix object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient array.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix.
 *
 *  @param a Object to be initialized.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_amatrix(pamatrix a, uint rows, uint cols);

/** @brief Initialize an @ref amatrix object to represent a submatrix.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref amatrix for the coefficient array, leading to a new
 *  matrix representing a submatrix of the source.
 *  Changes to the coefficients of the new matrix also change
 *  coefficients of the source matrix.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param a Object to be initialized.
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_sub_amatrix(pamatrix a, pamatrix src, uint rows, uint roff, uint cols,
    uint coff);

/** @brief Initialize an @ref amatrix object by a vector.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of an @ref avector in order to represent the matrix.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param a Object to be initialized.
 *  @param src Source vector.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_vec_amatrix(pamatrix a, pavector src, uint rows, uint cols);

/** @brief Initialize an @ref amatrix object using a given array for
 *  the coefficients.
 *
 *  Sets up the components of the object and uses the given array to
 *  represent the coefficients.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param a Object to be initialized.
 *  @param src Source array, should contain at least <tt>rows * cols</tt> elements.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_pointer_amatrix(pamatrix a, pfield src, uint rows, uint cols);

/** @brief Initialize an @ref amatrix object and set it to zero.
 *
 *  Sets up the components of the object, allocates storage for the
 *  coefficient array, and sets it to zero.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix.
 *
 *  @param a Object to be initialized.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_zero_amatrix(pamatrix a, uint rows, uint cols);

/** @brief Initialize an @ref amatrix object and set it to identity.
 *
 *  Sets up the components of the object, allocates storage for the
 *  coefficient array, sets the diagonal to one and all off-diagonal
 *  entries to zero.
 *  For square matrices, this yields an identity matrix.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix.
 *
 *  @param a Object to be initialized.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_identity_amatrix(pamatrix a, uint rows, uint cols);

/** @brief Uninitialize an @ref amatrix object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param a Object to be uninitialized. */
HEADER_PREFIX void
uninit_amatrix(pamatrix a);

/** @brief Create a new @ref amatrix object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_amatrix(uint rows, uint cols);

/** @brief Create a new @ref amatrix object representing a submatrix.
 *
 *  Allocates storage for the object, but uses part of the storage
 *  of another @ref amatrix object to keep the coefficients.
 *  Since the leading dimension of the source matrix is used, this
 *  allows us to work with a rectangular submatrix of the original
 *  matrix.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_sub_amatrix(pamatrix src, uint rows, uint roff, uint cols, uint coff);

/** @brief Create a new @ref amatrix object using a given array for
 *  the coefficients.
 *
 *  @param src Source array, should contain at least <tt>rows * cols</tt> elements.
 *  @param rows Number of rows for the new matrix.
 *  @param cols Number of columns for the new matrix.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_pointer_amatrix(field *src, uint rows, uint cols);

/** @brief Create a new @ref amatrix object representing a zero matrix.
 *
 *  Allocates storage for the object and sets all coefficients to
 *  zero.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_zero_amatrix(uint rows, uint cols);

/** @brief Create a new @ref amatrix object representing the identity.
 *
 *  Allocates storage for the object and sets diagonal entries to one
 *  and off-diagonal entries to zero.
 *  For square matrices, this yields an identity matrix.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_identity_amatrix(uint rows, uint cols);

/** @brief Delete an @ref amatrix object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they will otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param a Object to be deleted. */
HEADER_PREFIX void
del_amatrix(pamatrix a);

/** @brief Change the dimensions of an @ref amatrix object without
 *  preserving its coefficients.
 *
 *  Allocates new storage for the coefficients and releases the
 *  old storage.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param a Matrix to be resized.
 *  @param rows New number of rows.
 *  @param cols New number of columns. */
HEADER_PREFIX void
resize_amatrix(pamatrix a, uint rows, uint cols);

/** @brief Change the dimensions of an @ref amatrix object while
 *  preserving as many of its coefficients as possible.
 *
 *  Allocates new storage for the coefficients and copies as many
 *  of the old coefficients into it.
 *  If there are more rows or columns in the new matrix, the additional
 *  coefficients are left unintialized.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param a Matrix to be resized.
 *  @param rows New number of rows.
 *  @param cols New number of columns. */
HEADER_PREFIX void
resizecopy_amatrix(pamatrix a, uint rows, uint cols);

/* ------------------------------------------------------------
 Access methods
 ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX field
getentry_amatrix(pcamatrix, uint, uint) __attribute__ ((const,unused));
INLINE_PREFIX void
setentry_amatrix(pamatrix, uint, uint, field) __attribute__((unused));
INLINE_PREFIX field
addentry_amatrix(pamatrix, uint, uint, field) __attribute__((unused));
#endif

/** @brief Read a matrix entry @f$a_{ij}@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @returns Matrix entry @f$a_{ij}@f$. */
INLINE_PREFIX field getentry_amatrix(pcamatrix a, uint row, uint col) {
  longindex lda = a->ld;
#ifdef FULL_DEBUG
  assert(row < a->rows);
  assert(col < a->cols);
#endif

  return a->a[row + lda * col];
}

/** @brief Set a matrix entry, @f$a_{ij}\gets x@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x New value of @f$a_{ij}@f$. */
INLINE_PREFIX void setentry_amatrix(pamatrix a, uint row, uint col, field x) {
  longindex lda = a->ld;
#ifdef FULL_DEBUG
  assert(row < a->rows);
  assert(col < a->cols);
#endif

  a->a[row + lda * col] = x;
}

/** @brief Add to a matrix entry, @f$a_{ij} \gets a_{ij} + x@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x Summand.
 *  @returns New value of @f$a_{ij}@f$. */
INLINE_PREFIX field addentry_amatrix(pamatrix a, uint row, uint col, field x) {
  longindex lda = a->ld;
#ifdef FULL_DEBUG
  assert(row < a->rows);
  assert(col < a->cols);
#endif

  return (a->a[row + lda * col] += x);
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

/** @brief Get number of currently initialized @ref amatrix objects.
 *
 *  Calls to initialization functions like @ref init_amatrix and
 *  constructors like @ref new_amatrix increase an internal counter,
 *  while @ref uninit_amatrix and @ref del_amatrix decrease it.
 *
 *  @remark Use this function to check whether a program correctly cleans
 *  up temporary variables.
 *
 *  @returns Number of currently initialized @ref amatrix objects. */
HEADER_PREFIX uint
getactives_amatrix();

/** @brief Get size of a given @ref amatrix object.
 *
 *  Computes the size of the @ref amatrix object and the storage
 *  allocated for the coefficients.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_amatrix), no coefficient
 *  storage is added.
 *
 *  @param a Matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_amatrix(pcamatrix a);

/** @brief Get heap size of a given @ref amatrix object.
 *
 *  Computes the size of storage allocated for the coefficients
 *  on the heap, but not for the @ref amatrix object itself.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_amatrix), no storage
 *  is required.
 *
 *  @param a Matrix object.
 *  @returns Size of allocated heap storage in bytes. */
HEADER_PREFIX size_t
getsize_heap_amatrix(pcamatrix a);

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

/** @brief Set a matrix to zero.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
clear_amatrix(pamatrix a);

/** @brief Set the lower triangular part of a matrix to zero.
 *
 *  @param a Target matrix.
 *  @param strict If set, only the <em>strict</em> lower triangular
 *     part should be cleared. */
HEADER_PREFIX void
clear_lower_amatrix(pamatrix a, bool strict);

/** @brief Set the upper triangular part of a matrix to zero.
 *
 *  @param a Target matrix.
 *  @param strict If set, only the <em>strict</em> upper triangular
 *     part should be cleared. */
HEADER_PREFIX void
clear_upper_amatrix(pamatrix a, bool strict);

/** @brief Set a matrix to identity.
 *
 *  Sets all diagonal entries to one and all off-diagonal entries to zero.
 *  For square matrices, this yields an identity matrix.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
identity_amatrix(pamatrix a);

/** @brief Fill a matrix with random values.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
random_amatrix(pamatrix a);

/** @brief Fill a square matrix with random values and ensure that it is invertible
 *
 *  First the matrix is filled with random values, then the diagonal
 *  elements are set to @f$a_{ii} \gets \alpha + \sum_{j=1}^n |a_{ij}|@f$,
 *  ensuring that @f$A@f$ is diagonal dominant and therefore invertible.
 *
 *  @param a Target matrix.
 *  @param alpha Diagonal shift @f$\alpha@f$, should be a positive number. */
HEADER_PREFIX void
random_invertible_amatrix(pamatrix a, real alpha);

/** @brief Fill a matrix with random values and ensure that it is self-adjoint.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
random_selfadjoint_amatrix(pamatrix a);

/** @brief Fill a matrix with random values and ensure that it is positive definite.
 *
 *  First the matrix is filled with random values, ensuring that it becomes
 *  self-adjoint.
 *  Then the diagonal elements are set to
 *  @f$a_{ii} \gets \alpha + \sum_{j=1}^n |a_{ij}|@f$, ensuring that
 *  @f$A@f$ is diagonal dominant and therefore positiv definite.
 *
 *  @param a Target matrix.
 *  @param alpha Diagonal shift @f$\alpha@f$, should be a positiv number. */
HEADER_PREFIX void
random_spd_amatrix(pamatrix a, real alpha);

/** @brief Copy a matrix into another matrix, @f$B \gets A@f$ or
 *  @f$B \gets A^*@f$.
 *
 *  The numbers of rows and columns have to match.
 *
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Source matrix.
 *  @param b Target matrix. */
HEADER_PREFIX void
copy_amatrix(bool atrans, pcamatrix a, pamatrix b);

/** @brief Copy a matrix into another matrix, permuting the columns,
 *  i.e., @f$B \gets A P_\pi@f$ or @f$B \gets (A P_\pi)^*@f$.
 *
 *  The numbers of rows and columns have to match.
 *
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Source matrix.
 *  @param colpiv Array of dimension <tt>a->cols</tt> containing the
 *     column indices.
 *  @param b Target matrix. */
HEADER_PREFIX void
copy_colpiv_amatrix(bool atrans, pcamatrix a, const uint *colpiv,
		    pamatrix b);

/** @brief Create a duplicate of an existing @ref amatrix.
 *
 *  @param src Matrix to be duplicated.
 *  @returns Copy of <tt>src</tt>. */
HEADER_PREFIX pamatrix
clone_amatrix(pcamatrix src);

/** @brief Copy a matrix into another matrix, @f$B \gets A@f$ or
 *  @f$B \gets A^*@f$.
 *
 *  If @f$B@f$ is smaller than @f$A@f$ (or @f$A^*@f$), only the upper
 *  left part of @f$A@f$ is copied.
 *  If @f$A@f$ (or @f$A^*@f$) is smaller than @f$B@f$, only the upper
 *  left part of @f$B@f$ is filled.
 *
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Source matrix.
 *  @param b Target matrix. */
HEADER_PREFIX void
copy_sub_amatrix(bool atrans, pcamatrix a, pamatrix b);

/** @brief Print a matrix.
 *
 *  @param a Matrix object. */
HEADER_PREFIX void
print_amatrix(pcamatrix a);

/** @brief Print a matrix in Matlab format.
 *
 *  @param a Matrix object. */
HEADER_PREFIX void
print_matlab_amatrix(pcamatrix a);

/** @brief Check whether a matrix @f$A@f$ or its adjoint @f$A^*@f$ is isometric.
 *
 *  Compute either @f$I-A^*A@f$ or @f$I-AA^*@f$.
 *  If the return value is small, @f$A@f$ or @f$A^*@f$ are isometric
 *  matrices.
 *  For square matrices, this also means that they are orthogonal.
 *
 *  @param atrans Set if @f$A^*@f$ is to be checked, otherwise @f$A@f$ is used.
 *  @param a Matrix @f$A@f$.
 *  @returns @f$I-A^*A@f$ if <tt>atrans==false</tt> and @f$I-AA^*@f$ otherwise. */
HEADER_PREFIX real
check_ortho_amatrix(bool atrans, pcamatrix a);

/* ------------------------------------------------------------
 Basic linear algebra
 ------------------------------------------------------------ */

/** @brief Scale a matrix @f$A@f$ by a factor @f$\alpha@f$,
 *  @f$A \gets \alpha A@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Target matrix @f$A@f$. */
HEADER_PREFIX void
scale_amatrix(field alpha, pamatrix a);

/**
 * @brief compute the complex conjugate @f$ \bar A @f$ of a matrix @f$ A @f$.
 *
 * @param a Matrix that will be conjugated.
 */
HEADER_PREFIX void
conjugate_amatrix(pamatrix a);

/** @brief Compute the Frobenius inner product
 *  @f$\langle A, B \rangle_F@f$ of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The Frobenius inner product is given by
 *  @f$\langle A, B \rangle_F = \sum_{i,j} \bar a_{ij} b_{ij}@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @param b Matrix @f$B@f$.
 *  @returns Frobenius inner product @f$\langle A, B \rangle_F@f$. */
HEADER_PREFIX field
dotprod_amatrix(pcamatrix a, pcamatrix b);

/** @brief Approximate the spectral norm @f$\|A\|_2@f$ of a matrix @f$A@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$A^* A@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @returns Approximation of @f$\|A\|_2@f$. */
HEADER_PREFIX real
norm2_amatrix(pcamatrix A);

/** @brief Compute the Frobenius norm @f$\|A\|_F@f$ of a matrix @f$A@f$.
 *
 *  The Frobenius norm is given by
 *  @f$\|A\|_F = \left( \sum_{i,j} |a_{ij}|^2 \right)^{1/2}@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @returns Frobenius norm @f$\|A\|_F@f$. */
HEADER_PREFIX real
normfrob_amatrix(pcamatrix a);

/** @brief Compute the squared Frobenius norm @f$\|A\|_F^2@f$ of a matrix @f$A@f$.
 *
 *  The Frobenius norm is given by
 *  @f$\|A\|_F^2 = \sum_{i,j} |a_{ij}|^2@f$.
 *
 *  @param a Matrix @f$A@f$.
 *  @returns Squared Frobenius norm @f$\|A\|_F^2@f$. */
HEADER_PREFIX real
normfrob2_amatrix(pcamatrix a);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Dense matrix @f$A@f$.
 *  @param b Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix(pcamatrix a, pcamatrix b);

/** @brief Multiply a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix @f$A@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_amatrix_avector(field alpha, pcamatrix a, pcavector src, pavector trg);

/** @brief Multiply the adjoint of a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint @f$A^*@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_amatrix_avector(field alpha, pcamatrix a, pcavector src,
    pavector trg);

/** @brief Multiply a matrix @f$A@f$ or its adjoint @f$A^*@f$ by a
 *  vector, @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_amatrix_avector(field alpha, bool atrans, pcamatrix a, pcavector src,
    pavector trg);

/** @brief Add two matrices,
 *  @f$B \gets B + \alpha A@f$ or @f$B \gets B + \alpha A^*@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be added instead of @f$A@f$.
 *  @param a Source matrix @f$A@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix(field alpha, bool atrans, pcamatrix a, pamatrix b);

/** @brief Multiply two matrices,
 *  @f$C \gets C + \alpha A B@f$, @f$C \gets C + \alpha A^* B@f$,
 *  @f$C \gets C + \alpha A B^*@f$ or @f$C \gets C + \alpha A^* B^*@f$.
 *
 *  The matrices @f$A@f$ (or @f$A^*@f$) and @f$B@f$ (or @f$B^*@f$)
 *  are multiplied, the result is scaled by @f$\alpha@f$ and added
 *  to the target matrix @f$C@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Left factor @f$A@f$.
 *  @param btrans Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param b Right factor @f$B@f$.
 *  @param c Target matrix @f$C@f$. */
HEADER_PREFIX void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans, pcamatrix b,
    pamatrix c);

/** @brief Multiply a matrix by a bidiagonal matrix,
 *  @f$A \gets \alpha A L@f$ or @f$A \gets \alpha L^* A@f$.
 *
 *  The matrix @f$A@f$ (or @f$A^*@f$) is multiplied by the
 *  lower bidiagonal matrix @f$L@f$ described by the diagonal vector @f$d@f$
 *  and the subdiagonal vector @f$l@f$, the result is scaled by
 *  @f$\alpha@f$ and written back into @f$A@f$, i.e.,
 *  @f$a_{ij} \gets \alpha (a_{ij} d_j + a_{i,j+1} l_j)@f$ or
 *  @f$a_{ij} \gets \alpha (\bar d_i a_{ij} + \bar l_i a_{i,j+1})@f$
 *  (adding zeros to @f$b@f$ for the last column or row).
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param d Diagonal entries of @f$L@f$.
 *  @param l Subdiagonal entries of @f$L@f$. */
HEADER_PREFIX void
bidiagmul_amatrix(field alpha, bool atrans, pamatrix a, pcavector d,
    pcavector l);

/** @} */

#endif
