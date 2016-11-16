/* ------------------------------------------------------------
 * This is the file "rkmatrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

/** @file rkmatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef RKMATRIX_H
#define RKMATRIX_H

/** @defgroup rkmatrix rkmatrix
 *  @brief Representation of a low-rank matrix in factorized form
 *  @f$R = A B^*@f$.
 *
 *  The @ref rkmatrix class is used to approximate admissible blocks
 *  in hierarchical matrices.
 *  @{ */

/** @brief Representation of a low-rank matrix in factorized form. */
typedef struct _rkmatrix rkmatrix;

/** @brief Pointer to @ref rkmatrix object. */
typedef rkmatrix *prkmatrix;

/** @brief Pointer to constant @ref rkmatrix object. */
typedef const rkmatrix *pcrkmatrix;

#include "amatrix.h"
#include "krylov.h"
#include "settings.h"

/** @brief Representation of a low-rank matrix in factorized form
 *  @f$R = A B^*@f$. */
struct _rkmatrix {
  /** Row factor @f$A@f$. */
  amatrix A;
  /** Column factor @f$B@f$. */
  amatrix B;

  /** Maximal rank, i.e., number of columns of @f$A@f$ and @f$B@f$. */
  uint k;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Initialize an @ref rkmatrix object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the matrices @f$A@f$ and @f$B@f$.
 *
 *  @remark Should always be matched by a call to @ref uninit_rkmatrix.
 *
 *  @param r Object to be initialized.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @param k Rank.
 *  @returns Initialized @ref rkmatrix object. */
HEADER_PREFIX prkmatrix
init_rkmatrix(prkmatrix r, uint rows, uint cols, uint k);

/** @brief Initialize an @ref rkmatrix object to represent a submatrix.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref rkmatrix for the matrix factors, leading to a new
 *  matrix representing a submatrix of the source.
 *
 *  @remark The function returns a pointer to a constant object in order
 *  to discourage changing the submatrix. Due to its factorized
 *  representation, changing the submatrix would lead to non-local
 *  changes to the entire matrix, and this is only rarely a good idea.
 *
 *  @remark Should always be matched by a call to @ref uninit_rkmatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param r Object to be initialized.
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns Initialized @ref rkmatrix object. */
HEADER_PREFIX pcrkmatrix
init_sub_rkmatrix(prkmatrix r, pcrkmatrix src, uint rows, uint roff, uint cols,
    uint coff);

/** @brief Uninitialize an @ref rkmatrix object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param r Object to be uninitialized. */
HEADER_PREFIX void
uninit_rkmatrix(prkmatrix r);

/** @brief Create a new @ref rkmatrix object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_rkmatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @param k Rank.
 *  @returns New @ref rkmatrix object. */
HEADER_PREFIX prkmatrix
new_rkmatrix(uint rows, uint cols, uint k);

/** @brief Create a new @ref rkmatrix object representing a submatrix.
 *
 *  Allocates storage for the object and uses part of the storage
 *  of another @ref rkmatrix for the matrix factors, leading to a new
 *  matrix representing a submatrix of the source.
 *
 *  @remark The function returns a pointer to a constant object in order
 *  to discourage changing the submatrix. Due to its factorized
 *  representation, changing the submatrix would lead to non-local
 *  changes to the entire matrix, and this is only rarely a good idea.
 *
 *  @remark Should always be matched by a call to @ref del_rkmatrix.
 *
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns New @ref rkmatrix object. */
HEADER_PREFIX pcrkmatrix
new_sub_rkmatrix(pcrkmatrix src, uint rows, uint roff, uint cols, uint coff);

/** @brief Delete an @ref rkmatrix object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they will otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param r Object to be deleted. */
HEADER_PREFIX void
del_rkmatrix(prkmatrix r);

/* ------------------------------------------------------------
 * Change rank
 * ------------------------------------------------------------ */

/** @brief Change the rank of an @ref rkmatrix.
 *
 *  This function resizes the matrices @f$A@f$ and @f$B@f$
 *  @e without preserving their contents.
 *  The caller is responsible for ensuring that these matrices
 *  are initialized properly after they have been resized.
 *
 *  @param r Target matrix.
 *  @param k New rank. */
void
setrank_rkmatrix(prkmatrix r, uint k);

/** @brief Change the size of an @ref rkmatrix.
 *
 *  This function resizes the matrices @f$A@f$ and @f$B@f$
 *  @e without preserving their contents.
 *  The caller is responsible for ensuring that these matrices
 *  are initialized properly after they have been resized.
 *
 *  @param r Target matrix.
 *  @param rows New number of rows.
 *  @param cols New number of columns.
 *  @param k New rank. */
void
resize_rkmatrix(prkmatrix r, uint rows, uint cols, uint k);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref rkmatrix object.
 *
 *  Computes the size of the @ref rkmatrix object and the storage
 *  allocated for the coefficients.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref init_sub_rkmatrix), no coefficient
 *  storage is added.
 *
 *  @param r Matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_rkmatrix(pcrkmatrix r);

/** @brief Get heap size of a given @ref rkmatrix object.
 *
 *  Computes the size of storage allocated for the coefficients
 *  on the heap, but not for the @ref rkmatrix object itself.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_rkmatrix), no storage
 *  is required.
 *
 *  @param r Matrix object.
 *  @returns Size of allocated heap storage in bytes. */
HEADER_PREFIX size_t
getsize_heap_rkmatrix(pcrkmatrix r);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Create a copy of an @ref rkmatrix.
 *
 *  Creates a new @ref rkmatrix with identical parameters and
 *  copies the contents of @f$A@f$ and @f$B@f$.
 *
 *  @param r Source @ref rkmatrix.
 *  @returns Full copy of <tt>r</tt>. */
HEADER_PREFIX prkmatrix
clone_rkmatrix(pcrkmatrix r);

/** @brief Copy an @ref rkmatrix into another @ref rkmatrix.
 *
 *  Copies the contents of <tt>a->A</tt> and <tt>a->B</tt> to
 *  <tt>b->A</tt> and <tt>b->B</tt>.
 *
 *  @param atrans Determines wether @f$A@f$ or @f$A^*@f$ should be copied.
 *  @param a Source @ref rkmatrix.
 *  @param b Target @ref rkmatrix.
 */
HEADER_PREFIX void
copy_rkmatrix(bool atrans, pcrkmatrix a, prkmatrix b);

/** @brief Scale an @ref rkmatrix by a factor.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param r Target @ref rkmatrix @f$R@f$, will be overwritten
 *           by @f$\alpha R@f$. */
HEADER_PREFIX void
scale_rkmatrix(field alpha, prkmatrix r);

/** @brief Fill an @ref rkmatrix with random values.
 *
 *  @param r Target matrix.
 *  @param kmax Maximal rank. */
HEADER_PREFIX void
random_rkmatrix(prkmatrix r, uint kmax);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Multiply a low-rank matrix @f$R@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha R x@f$.
 *
 *  The matrix @f$R@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param r Matrix @f$R@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_rkmatrix_avector(field alpha, pcrkmatrix r, pcavector x, pavector y);

/** @brief Multiply the adjoint of a low-rank matrix @f$R@f$ by a
 *  vector @f$x@f$, @f$y \gets y + \alpha R^* x@f$.
 *
 *  The adjoint @f$R^*@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param r Matrix @f$R@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_rkmatrix_avector(field alpha, pcrkmatrix r, pcavector x,
    pavector y);

/** @brief Multiply a low-rank matrix @f$R@f$ or its adjoint @f$R^*@f$ by
 *  a vector, @f$y \gets y + \alpha R x@f$ or @f$y \gets y + \alpha R^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param rtrans Set if @f$R^*@f$ is to be used instead of @f$R@f$.
 *  @param r Matrix @f$R@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_rkmatrix_avector(field alpha, bool rtrans, pcrkmatrix r, pcavector x,
    pavector y);

/* ------------------------------------------------------------
 * Spectral norm
 * ------------------------------------------------------------ */

/** @brief Approximate the spectral norm @f$\|R\|_2@f$ of a matrix @f$R@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$R^* R@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param R Low rank matrix @f$R@f$.
 *  @returns Approximation of @f$\|R\|_2@f$. */
HEADER_PREFIX real
norm2_rkmatrix(pcrkmatrix R);

/** @brief Approximate the spectral norm @f$\|R_A-R_B\|_2@f$ of the difference
 *  of two low-rank matrices @f$R_A@f$ and @f$R_B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(R_A-R_B)^* (R_A-R_B)@f$ and computing the
 *  square root of the resulting eigenvalue approximation.
 *
 *  @param a Low rank matrix @f$R_A@f$.
 *  @param b Low rank matrix @f$R_B@f$.
 *  @returns Approximation of @f$\|R_A-R_B\|_2@f$. */
HEADER_PREFIX real
norm2diff_rkmatrix(pcrkmatrix a, pcrkmatrix b);

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

#ifndef RKMATRIX_COMPLETE
#define RKMATRIX_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_rkmatrix(pcrkmatrix) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_rkmatrix(pcrkmatrix) __attribute__ ((const,unused));
INLINE_PREFIX uint
getrank_rkmatrix(pcrkmatrix r) __attribute__ ((const, unused));
INLINE_PREFIX pamatrix
getA_rkmatrix(prkmatrix r) __attribute__ ((unused));
INLINE_PREFIX pamatrix
getB_rkmatrix(prkmatrix r) __attribute__ ((unused));
#endif

/** @brief Get the number of rows of an @ref rkmatrix @f$R=A B^*@f$.
 *
 *  @param r Matrix @f$R@f$.
 *  @return Number of rows of @f$R@f$, i.e., number of rows of @f$A@f$. */
INLINE_PREFIX uint
getrows_rkmatrix(pcrkmatrix r) {
  return r->A.rows;
}

/** @brief Get the number of columns of an @ref rkmatrix @f$R=A B^*@f$.
 *
 *  @param r Matrix @f$R@f$.
 *  @returns Number of columns of @f$R@f$, i.e., number of rows of @f$B@f$. */
INLINE_PREFIX uint
getcols_rkmatrix(pcrkmatrix r) {
  return r->B.rows;
}

/** @brief Get the rank of an @ref rkmatrix @f$R=A B^*@f$.
 *
 *  @param r Matrix @f$R@f$.
 *  @returns Maximal rank of @f$R@f$, i.e., number of columns of @f$A@f$ and @f$B@f$. */
INLINE_PREFIX uint
getrank_rkmatrix(pcrkmatrix r) {
  assert(r->A.cols == r->B.cols);

  return r->A.cols;
}

/** @brief Get the factor A of an @ref rkmatrix @f$R=A B^*@f$.
 *
 * @param r Matrix @f$R@f$.
 * @returns Factor @f$A@f$. */
INLINE_PREFIX pamatrix
getA_rkmatrix(prkmatrix r) {
  return &r->A;
}

/** @brief Get the factor B of an @ref rkmatrix @f$R=A B^*@f$.
 *
 * @param r Matrix @f$R@f$.
 * @returns Factor @f$A@f$. */
INLINE_PREFIX pamatrix
getB_rkmatrix(prkmatrix r) {
  return &r->B;
}

#endif

/** @} */
