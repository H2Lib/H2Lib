
/* ------------------------------------------------------------
   This is the file "hmatrix.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

/** @file hmatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef HMATRIX_H
#define HMATRIX_H

/** @defgroup hmatrix hmatrix
 *  @brief Representation of a hierarchical matrix.
 *
 *  The @ref hmatrix class is used to represent hierarchical matrices
 *  with arbitrary block structures and arbitrary rank distributions.
 *  @{ */

/** @brief Representation of a hierarchical matrix. */
typedef struct _hmatrix hmatrix;

/** @brief Pointer to a @ref hmatrix object. */
typedef hmatrix *phmatrix;

/** @brief Pointer to constant @ref hmatrix object. */
typedef const hmatrix *pchmatrix;

#ifdef USE_CAIRO
#include <cairo.h>
#endif

#include "amatrix.h"
#include "factorizations.h"
#include "block.h"
#include "rkmatrix.h"
#include "settings.h"
#include "eigensolvers.h"

/** @brief Representation of @f$\mathcal{H}@f$-matrices.
 *
 *  @f$\mathcal{H}@f$-matrices are represented recursively:
 *  an @ref hmatrix object can be either an @ref amatrix,
 *  an @ref rkmatrix or divided into submatrices represented
 *  again by @ref hmatrix objects. */
struct _hmatrix {
  /** @brief Row cluster */
  pccluster rc;
  /** @brief Column cluster */
  pccluster cc;

  /** @brief Low-rank matrix, for admissible leaves. */
  prkmatrix r;

  /** @brief Standard matrix, for inadmissible leaves. */
  pamatrix f;
  
  /** @brief Submatrices. */
  phmatrix *son;
  /** @brief Number of block rows. */
  uint rsons;
  /** @brief Number of block columns. */
  uint csons;

  /** @brief Number of references to this @ref hmatrix. */
  uint refs;
  /** @brief Number of descendants in matrix tree. */
  uint desc;
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

HEADER_PREFIX phmatrix
init_hmatrix(phmatrix hm, pccluster rc, pccluster cc);

HEADER_PREFIX void
uninit_hmatrix(phmatrix hm);

/** @brief Create a new @ref hmatrix object.
 *
 *  Allocates storage for the object and sets the matrix pointers
 *  to NULL, representing a zero matrix.
 *
 *  @remark Should always be matched by a call to @ref del_hmatrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_hmatrix) is deleted.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @returns New @ref hmatrix object. */
HEADER_PREFIX phmatrix
new_hmatrix(pccluster rc, pccluster cc);

/** @brief Create a new @ref hmatrix object representing a
 *  low-rank matrix.
 *
 *  Allocates storage for the object containing an @ref rkmatrix
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_hmatrix.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param k Rank.
 *  @returns New @ref hmatrix object containing a new @ref rkmatrix. */
HEADER_PREFIX phmatrix
new_rk_hmatrix(pccluster rc, pccluster cc, uint k);

/** @brief Create a new @ref hmatrix object representing a
 *  standard dense matrix.
 *
 *  Allocates storage for the object containing an @ref amatrix
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_hmatrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_hmatrix) is deleted.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @returns New @ref hmatrix object containing a new @ref amatrix. */
HEADER_PREFIX phmatrix
new_full_hmatrix(pccluster rc, pccluster cc);

/** @brief Create a new @ref hmatrix object representing a
 *  subdivided matrix.
 *
 *  Allocates storage for the object representing a matrix with
 *  submatrices. The submatrices are initialized with NULL pointers
 *  and have to be set up with @ref ref_hmatrix.
 *
 *  @remark Should always be matched by a call to @ref del_hmatrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_hmatrix) is deleted.
 *
 *  @param rc Row cluster.
 *  @param cc Column cluster.
 *  @param rsons Number of block rows.
 *  @param csons Number of block columns.
 *  @returns New @ref hmatrix object with submatrices. */
HEADER_PREFIX phmatrix
new_super_hmatrix(pccluster rc, pccluster cc,
		  uint rsons, uint csons);

/** @brief Creates a clone of an existing @ref _hmatrix "hmatrix".
 *
 *  This function creates a clone of an existing @ref _hmatrix "hmatrix".
 *  I.e. a new @ref _hmatrix "hmatrix" object is created with the same block
 *  structure as the input matrix. Further all @ref _rkmatrix "rank k matrices"
 *  and @ref _amatrix "dense matrices" are copied to the clone aswell.
 *
 *  @param src @ref _hmatrix "Hmatrix" to be cloned.
 *  @returns A clone of <tt>src</tt> is returned. */
HEADER_PREFIX phmatrix
clone_hmatrix(pchmatrix src);

/** @brief Clones the structure of an existing @ref _hmatrix "hmatrix".
 *
 *  This function clones the structure of an existing @ref _hmatrix "hmatrix".
 *  I.e. a new @ref _hmatrix "hmatrix" object is created with the same block
 *  structure and local ranks as the input matrix.
 *  But the matrix entries from <tt>src</tt> are not copied out to the
 *  resulting matrix.
 *
 *  @param src @ref _hmatrix "Hmatrix" to be cloned.
 *  @returns A clone of <tt>src</tt> is returned. */
HEADER_PREFIX phmatrix
clonestructure_hmatrix(pchmatrix src);

/** @brief Complete the initialisation of a @ref hmatrix object.
 *
 *  Complete the initialisation of the @ref hmatrix object after all sons
 *  have been initialised. 
 *  The number of the descendants of the H-matrix is computed from
 *  the <tt>desc</tt> fields of its sons.
 * 
 *  @param hm H-matrix to be completed. */
HEADER_PREFIX void
update_hmatrix(phmatrix hm);

/** @brief Delete an @ref hmatrix object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref hmatrix contains pointers to submatrices,
 *  the submatrices are released by @ref unref_hmatrix.
 *
 *  Only objects with <tt>hm->refs==0</tt> may be deleted.
 *
 *  @param hm Object to be deleted. */
HEADER_PREFIX void
del_hmatrix(phmatrix hm);

/* ------------------------------------------------------------
   Reference counting
   ------------------------------------------------------------ */

/** @brief Set a pointer to an @ref hmatrix object, increase its
 *  reference counter, and decrease reference counter of original
 *  pointer target.
 *
 *  @param ptr Pointer to the @ref phmatrix variable that will be changed.
 *  @param hm @ref hmatrix that will be referenced. */
HEADER_PREFIX void
ref_hmatrix(phmatrix *ptr, phmatrix hm);


/** @brief Reduce the reference counter of an @ref hmatrix object.
 *
 *  If the reference counter reaches zero, the object is deleted.
 *
 *  @remark Use @ref ref_hmatrix with <tt>hm=NULL</tt> instead, since this
 *  guarantees that the pointer is properly invalidated.
 *
 *  @param hm @ref hmatrix that will be unreferenced. */
HEADER_PREFIX void
unref_hmatrix(phmatrix hm);

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

/** @brief Get size of a given @ref hmatrix object.
 *
 *  @param hm @f$\mathcal{H}@f$-matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_hmatrix(pchmatrix hm);

/** @brief Get size of the nearfield part of a given @ref hmatrix object.
 *
 *  @param hm @f$\mathcal{H}@f$-matrix object.
 *  @returns Size of allocated storage for nearfield in bytes. */
HEADER_PREFIX size_t
getnearsize_hmatrix(pchmatrix hm);

/** @brief Get size of the farfield part of a given @ref hmatrix object.
 *
 *  @param hm @f$\mathcal{H}@f$-matrix object.
 *  @returns Size of allocated storage for farfield in bytes. */
HEADER_PREFIX size_t
getfarsize_hmatrix(pchmatrix hm);

/* ------------------------------------------------------------
   Simple utility functions
   ------------------------------------------------------------ */

/** @brief Set an @ref hmatrix to zero by clearing all far- and nearfield
 *  matrices.
 *
 *  @param hm Target matrix. */
HEADER_PREFIX void
clear_hmatrix(phmatrix hm);

/**
 * @brief Copy a matrix <tt>src</tt> to an existing matrix <tt>dst</tt>.
 *
 * It is assumed, that <tt>dst</tt> has exactly the same block structure
 * as <tt>src</tt>.
 *
 * @param src Source matrix.
 * @param trg Target matrix.
 */
HEADER_PREFIX void
copy_hmatrix(pchmatrix src, phmatrix trg);

/* ------------------------------------------------------------
   Build H-matrix based on block tree
   ------------------------------------------------------------ */

/** @brief Build an @ref hmatrix object from a @ref block tree using
 *  a given local rank.
 *
 *  @remark Submatrices for far- and nearfield leaves are created,
 *  but their coefficients are not initialized.
 *
 *  @param b Block tree.
 *  @param k Local rank.
 *  @returns New @ref hmatrix object. */
HEADER_PREFIX phmatrix
build_from_block_hmatrix(pcblock b, uint k);

/* ------------------------------------------------------------
   Matrix-vector multiplication
   ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_hmatrix_avector(field alpha, bool atrans, pchmatrix a,
	    pcavector x, pavector y);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param xp Source vector @f$x@f$ in cluster numbering
 *            with respect to <tt>hm->cc</tt>.
 *  @param yp Target vector @f$y@f$ in cluster numbering
 *            with respect to <tt>hm->rc</tt>. */
HEADER_PREFIX void
fastaddeval_hmatrix_avector(field alpha, pchmatrix hm, pcavector xp, pavector yp);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_hmatrix_avector(field alpha, pchmatrix hm,
		pcavector x, pavector y);

/** @brief Adjoint matrix-vector multiplication
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param xp Source vector @f$x@f$ in cluster numbering
 *            with respect to <tt>hm->rc</tt>.
 *  @param yp Target vector @f$y@f$ in cluster numbering
 *            with respect to <tt>hm->cc</tt>. */
HEADER_PREFIX void
fastaddevaltrans_hmatrix_avector(field alpha, pchmatrix hm,
			 pcavector xp, pavector yp);

/** @brief Adjoint matrix-vector multiplication
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_hmatrix_avector(field alpha, pchmatrix hm,
		     pcavector x, pavector y);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$ with symmetric matrix @f$A@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *  Only the lower triangular part of @f$A@f$ is used.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param xp Source vector @f$x@f$ in cluster numbering
 *            with respect to <tt>hm->cc</tt>.
 *  @param yp Target vector @f$y@f$ in cluster numbering
 *            with respect to <tt>hm->rc</tt>. */
HEADER_PREFIX void
fastaddevalsymm_hmatrix_avector(field alpha, pchmatrix hm,
				pcavector xp, pavector yp);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$ with symmetric matrix @f$A@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *  Only the lower triangular part of @f$A@f$ is used.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param hm Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevalsymm_hmatrix_avector(field alpha, pchmatrix hm,
			    pcavector x, pavector y);

/* ------------------------------------------------------------
   Enumeration by block number
   ------------------------------------------------------------ */

/** @brief Enumerate @ref hmatrix according to block tree.
 *
 *  The @ref hmatrix submatrices are enumerated in an array of
 *  size <tt>b->desc</tt>. The enumeration starts with <tt>0</tt> assigned to
 *  the root and then proceeds column-wise starting with <tt>b->sons[0]</tt>
 *  corresponding to the entries <tt>1</tt> to <tt>b->sons[0]->desc</tt>
 *  in the array and ending with <tt>b->sons[rsons*csons-1]</tt> corresponding
 *  to the last <tt>b->sons[rsons*csons-1]->desc</tt> entries.
 *
 *  @param b Block tree.
 *  @param hm Matrix matching the block structure given by <tt>b</tt>.
 *  @returns Array of size <tt>b->desc</tt> containing pointers to the
 *         @ref hmatrix objects corresponding to descendants of <tt>hm</tt>. */
HEADER_PREFIX phmatrix *
enumerate_hmatrix(pcblock b, phmatrix hm);

/* ------------------------------------------------------------
   Spectral norm
   ------------------------------------------------------------ */

/** @brief Approximate the spectral norm @f$\|A\|_2@f$ of a matrix @f$A@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$A^* A@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param a Matrix @f$A@f$.
 *  @returns Approximation of @f$\|A\|_2@f$. */
HEADER_PREFIX real
norm2_hmatrix(pchmatrix a);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a Matrix @f$A@f$.
 *  @param b Matrix @f$B@f$.
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
 *  @param a Matrix @f$A@f$.
 *  @param b Matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_hmatrix(pchmatrix a, pchmatrix b);

/* ------------------------------------------------------------
   File I/O
   ------------------------------------------------------------ */

/** @brief Write a matrix into an ASCII file in the old HLib format.
 *
 *  @param G Hierarchical matrix.
 *  @param filename Name of the target file. */
void
write_hlib_hmatrix(pchmatrix G, const char *filename);

/** @brief Read a matrix from an ASCII file in the old HLib format.
 *
 *  @remark Since HLib does not store cluster trees, this function
 *  has to reconstruct them from available data. It will always use
 *  one-dimensional trees with no permutation.
 *
 *  @param filename Name of the source file.
 *  @returns Hierarchical matrix reconstructed from file. */
phmatrix
read_hlib_hmatrix(const char *filename);

/* ------------------------------------------------------------
   Drawing
   ------------------------------------------------------------ */

#ifdef USE_CAIRO
HEADER_PREFIX void
draw_cairo_hmatrix(cairo_t *cr, pchmatrix hm,
		  bool storage, uint levels);
#endif

/** @} */

#endif
