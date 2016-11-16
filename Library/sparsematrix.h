
/* ------------------------------------------------------------
 * This is the file "sparsematrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2012
 * ------------------------------------------------------------ */

/** @file sparsematrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

/** @defgroup sparsematrix sparsematrix
 *  @brief Representation of a sparse matrix in compressed row format.
 *
 *  The @ref sparsematrix class is used to store sparse matrices
 *  efficiently by avoiding to store zero entries.
 *  The internal representation is the compressed row format:
 *  the @f$i@f$-th row is described by the entries @c row[i] to
 *  @c row[i+1]-1 in the arrays @c col and @c coeff.
 *  For each @f$a_{ij}\neq 0@f$ there is exactly one @f$k@f$ such that
 *  <tt>row[i]<=k<row[i+1]</tt>, @c col[k]=j and
 *  <tt>coeff[k]</tt>@f$=a_{ij}@f$.
 *  @{ */

/** @brief Representation of a sparse matrix in compressed row format. */
typedef struct _sparsematrix sparsematrix;

/** @brief Pointer to @ref sparsematrix object. */
typedef sparsematrix *psparsematrix;

/** @brief Pointer to constant @ref sparsematrix object. */
typedef const sparsematrix *pcsparsematrix;

#include "basic.h"
#include "avector.h"
#include "sparsepattern.h"
#include "krylov.h"

/** @brief Representation of a sparse matrix in compressed row format. */
struct _sparsematrix {
  /** @brief Number of rows. */
  uint rows;
  /** @brief Number of columns. */
  uint cols;
  /** @brief Number of non-zero entries. */
  uint nz;

  /** @brief Starting indices for row representations in @c col and
   *  @c coeff. */
  uint *row;
  /** @brief Column indices of non-zero entries. */
  uint *col;
  /** @brief Coefficients of non-zero entries. */
  pfield coeff;
};

/* ------------------------------------------------------------ *
 * Constructors and destructors                                 *
 * ------------------------------------------------------------ */

/** @brief Create a sparsematrix object without initializing its arrays.
 *
 *  @remark Should always be matched by a call to @ref del_sparsematrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @param nz Number of non-zero entries.
 *  @returns Allocated @ref sparsematrix object, arrays @c row,
 *     @c col and @c coeff are uninitialized. */
HEADER_PREFIX psparsematrix
new_raw_sparsematrix(uint rows, uint cols, uint nz);

/**
 * @brief Creates a new @ref _sparsematrix "sparsematrix" object and
 * initializes it to the idenity matrix @f$I \in \mathbb R^{n \times m}@f$.
 *
 * @param rows Number of rows is equal to @f$n@f$
 * @param cols Number of columns is equal to @f$m@f$
 * @return New identity sparsematrix
 */
HEADER_PREFIX psparsematrix
new_identity_sparsematrix(uint rows, uint cols);

/** @brief Create a sparsematrix based on a sparsepattern.
 *
 *  @remark Should always be matched by a call to @ref del_sparsematrix.
 *
 *  @param sp Description of matrix graph in a @ref sparsepattern object.
 *  @returns Fully initialized @ref sparsematrix object with given
 *     sparsity pattern and zero coefficients. */
HEADER_PREFIX psparsematrix
new_zero_sparsematrix(psparsepattern sp);

/** @brief Delete a @ref sparsematrix object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @param a Object to be deleted. */
HEADER_PREFIX void
del_sparsematrix(psparsematrix a);

/* ------------------------------------------------------------ *
 * Access methods
 * ------------------------------------------------------------ */

/** @brief Add to a matrix entry, @f$a_{ij} \gets a_{ij} + x@f$.
 *
 *  Only entries appearing in the sparsity pattern of the matrix
 *  are allowed.
 *
 *  @param a Target matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x Summand.
 *  @returns New value of @f$a_{ij}@f$. */
HEADER_PREFIX field
addentry_sparsematrix(psparsematrix a, uint row, uint col, field x);

/** @brief Set a matrix entry, @f$a_{ij} \gets x@f$.
 *
 *  Only entries appearing in the sparsity pattern of the matirx
 *  are allowed.
 *
 *  @param a Target matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x New value of @f$a_{ij}@f$. */
HEADER_PREFIX void
setentry_sparsematrix(psparsematrix a, uint row, uint col, field x);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref sparsematrix object.
 *
 *  Computes the size of the @ref sparsematrix object and the storage
 *  allocated for the coefficients.
 *
 *  @param a Target matrix.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_sparsematrix(pcsparsematrix a);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Sort non-zero entries to ensure that diagonal entries
 *     come first.
 *
 *  Reorder the entries of each row in @c col and @c coeff to
 *  place the diagonal entry first.
 *  Since many iterative solvers require us to handle this entry
 *  differently from all others, this optimization can improve
 *  the performance.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
sort_sparsematrix(psparsematrix a);

/** @brief Set a matrix to zero.
 *
 *  @param a Target matrix. */
HEADER_PREFIX void
clear_sparsematrix(psparsematrix a);

/** @brief Print a sparse matrix.
 *
 *  @param a Source matrix. */
HEADER_PREFIX void
print_sparsematrix(pcsparsematrix a);

/** @brief Print matrix to a Postscript file.
 *
 *  @param a Source matrix.
 *  @param filename Name of the target file.
 *  @param offset Offset added to all coordinates, e.g., to avoid
 *     boundary clipping with some printers. */
HEADER_PREFIX void
print_eps_sparsematrix(pcsparsematrix a, const char *filename, uint offset);

/* ------------------------------------------------------------
 * Basic linear algebra
 * ------------------------------------------------------------ */

/** @brief Multiply a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix @f$A@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_sparsematrix_avector(field alpha, pcsparsematrix a, pcavector x,
    pavector y);

/** @brief Multiply the adjoint of a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint @f$A^*@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_sparsematrix_avector(field alpha, pcsparsematrix a, pcavector x,
    pavector y);

/** @brief Multiply a matrix @f$A@f$ or its adjoint @f$A^*@f$ by a
 *  vector, @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param trans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_sparsematrix_avector(field alpha, bool trans, pcsparsematrix a, pcavector x,
    pavector y);

/** @brief Approximate the spectral norm @f$\|S\|_2@f$ of a matrix @f$S@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$S^* S@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param S Sparse matrix @f$S@f$.
 *  @returns Approximation of @f$\|S\|_2@f$. */
HEADER_PREFIX real
norm2_sparsematrix(pcsparsematrix S);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Sparse matrix @f$A@f$.
 *  @param b Sparse matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_sparsematrix(pcsparsematrix a, pcsparsematrix b);

/** @brief Add a @ref sparsematrix to an @ref amatrix,
 *  @f$B \gets B + \alpha A@f$ or @f$B \gets B + \alpha^* A@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be added instead of @f$A@f$.
 *  @param a Source matrix @f$A@f$.
 *  @param b Target matrix @f$B@f$. */
HEADER_PREFIX void
add_sparsematrix_amatrix(field alpha, bool atrans, pcsparsematrix a, pamatrix b);

/** @} */

#endif
