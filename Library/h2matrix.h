
/* ------------------------------------------------------------
 * This is the file "h2matrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file h2matrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef H2MATRIX_H
#define H2MATRIX_H

/** @defgroup h2matrix h2matrix
 *  @brief Representation of an @f$\mathcal{H}^2@f$-matrix.
 *
 *  The @ref h2matrix class is used to represent @f$\mathcal{H}^2@f$-matrices
 *  with arbitrary block structures and arbitrary cluster bases.
 *  @{ */

/** @brief Representation of an @f$\mathcal{H}^2@f$-matrix. */
typedef struct _h2matrix h2matrix;

/** @brief Pointer to @ref h2matrix object. */
typedef h2matrix *ph2matrix;

/** @brief Pointer to constant @ref h2matrix object. */
typedef const h2matrix *pch2matrix;

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "amatrix.h"
#include "krylov.h"
#include "block.h"
#include "hmatrix.h"
#include "uniform.h"
#include "clusterbasis.h"
#include "settings.h"

/** @brief Representation of @f$\mathcal{H}^2@f$-matrices.
 *
 *  @f$\mathcal{H}^2@f$-matrices are represented recursively:
 *  an @ref h2matrix object can be either an @ref amatrix,
 *  a @ref uniform matrix or divided into submatrices represented
 *  again by @ref h2matrix objects. */
struct _h2matrix {
  /** @brief Row cluster basis. */
  pclusterbasis rb;
  /** @brief Column cluster basis. */
  pclusterbasis cb;

  /** @brief Uniform matrix, for admissible leaves. */
  puniform u;

  /** @brief Standard matrix, for inadmissible leaves. */
  pamatrix f;

  /** @brief Submatrices. */
  ph2matrix *son;
  /** @brief Number of block rows. */
  uint rsons;
  /** @brief Number of block columns. */
  uint csons;

  /** @brief Number of references to this @ref h2matrix. */
  uint refs;
  /** @brief Number of descendants in matrix tree. */
  uint desc;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new @ref h2matrix object.
 *
 *  Allocates storage for the object and sets the matrix pointers
 *  to <tt>NULL,</tt> representing a zero matrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
new_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref h2matrix object representing a
 *  @ref uniform matrix.
 *
 *  Allocates storage for the object containing a @ref uniform
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object containing a new @ref uniform matrix. */
HEADER_PREFIX ph2matrix
new_uniform_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref h2matrix object representing a
 *  standard dense matrix.
 *
 *  Allocates storage for the object containing an @ref amatrix
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object containing a new @ref amatrix. */
HEADER_PREFIX ph2matrix
new_full_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref hmatrix object representing a
 *  subdivided matrix.
 *
 *  Allocates storage for the object representing a matrix with
 *  submatrices. The submatrices are initialized with NULL pointers
 *  and have to be set up with @ref ref_hmatrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @param rsons Number of block rows.
 *  @param csons Number of block columns.
 *  @returns New @ref h2matrix object with submatrices. */
HEADER_PREFIX ph2matrix
new_super_h2matrix(pclusterbasis rb, pclusterbasis cb, uint rsons, uint csons);

/** @brief Create a new @ref h2matrix object representing a zero matrix.
 *
 *  Allocates storage for the object representing a zero matrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
new_zero_h2matrix(pclusterbasis rb, pclusterbasis cb);

/**
 * @brief Builds a new @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the block structure of <tt>h2</tt>.
 *
 * @param h2 Input @ref h2matrix whose structure has to be cloned.
 * @param rb Row @ref clusterbasis for the clone.
 * @param cb Column @ref clusterbasis for the clone.
 * @return New @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the same block structure as <tt>h2</tt>.
 */
HEADER_PREFIX ph2matrix
clonestructure_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb);

/* clones a h2matrix, should be in the module h2matrix */
/**
 * @brief Builds a new @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the block structure of <tt>h2</tt>. Copies all coupling
 * matrices @f$ S_b,\ b \in \mathcal L_{\mathcal I \times \mathcal J}^+@f$ to
 * the clone.
 *
 * @param h2 Input @ref h2matrix which has to be cloned.
 * @param rb Row @ref clusterbasis for the clone.
 * @param cb Column @ref clusterbasis for the clone.
 * @return New @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the same data as <tt>h2</tt>.
 */
HEADER_PREFIX ph2matrix
clone_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb);

/** @brief Complete the initialisation of a @ref h2matrix object.
 *
 * Complete the initialisation of the @ref h2matrix object after all sons
 * have been initialised. 
 * The number of the descendants of the @f$\mathcal{H}^2@f$-matrix is
 * computed from the <tt>desc</tt> fields of its sons.
 * 
 * @param h2 @f$\mathcal{H}^2@f$-matrix to be completed. */
HEADER_PREFIX void
update_h2matrix(ph2matrix h2);

/** @brief Delete an @ref h2matrix object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref h2matrix contains pointers to submatrices,
 *  the submatrices are released by @ref unref_h2matrix.
 *
 *  Only objects with <tt>h2->refs==0</tt> may be deleted.
 *
 *  @param h2 Object to be deleted. */
HEADER_PREFIX void
del_h2matrix(ph2matrix h2);

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

/** @brief Set a pointer to an @ref h2matrix object, increase its
 *  reference counter, and decrease reference counter of original
 *  pointer target.
 *
 *  @param ptr Pointer to the @ref ph2matrix variable that will be changed.
 *  @param h2 @ref h2matrix that will be referenced. */
HEADER_PREFIX void
ref_h2matrix(ph2matrix *ptr, ph2matrix h2);

/** @brief Reduce the reference counter of an @ref h2matrix object.
 *
 *  If the reference counter reaches zero, the object is deleted.
 *
 *  @remark Use @ref ref_h2matrix with <tt>h2=NULL</tt> instead, since this
 *  guarantees that the pointer is properly invalidated.
 *
 *  @param h2 @ref h2matrix that will be unreferenced. */
HEADER_PREFIX void
unref_h2matrix(ph2matrix h2);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_h2matrix(pch2matrix h2);

/** @brief Get total size of a given @ref h2matrix object, including
 *  cluster bases.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
gettotalsize_h2matrix(pch2matrix h2);

/** @brief Get size of the nearfield part of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage for nearfield in bytes. */
HEADER_PREFIX size_t
getnearsize_h2matrix(pch2matrix h2);

/** @brief Get size of the farfield part of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage for farfield in bytes. */
HEADER_PREFIX size_t
getfarsize_h2matrix(pch2matrix h2);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Set an @ref h2matrix to zero by clearing all far- and nearfield
 *  matrices.
 *
 *  @param h2 Target matrix. */
HEADER_PREFIX void
clear_h2matrix(ph2matrix h2);

/** @brief Scale an @ref h2matrix by a factor.
 *
 * @param alpha Scaling factor @f$\alpha@f$.
 * @param h2 Target matrix @f$G@f$, will be overwritten by @f$\alpha G@f$. */
void
scale_h2matrix(field alpha, ph2matrix h2);

/** @brief Fill an @ref h2matrix with random coefficients.
 *
 * @param h2 Target matrix @f$G@f$. */
void
random_h2matrix(ph2matrix h2);

/* ------------------------------------------------------------
 * Build H^2-matrix based on block tree
 * ------------------------------------------------------------ */

/** @brief Build an @ref h2matrix object from a @ref block tree using
 *  given cluster bases.
 *
 *  @remark Submatrices for far- and nearfield leaves are created,
 *  but their coefficients are not initialized.
 *
 *  @param b Block tree.
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
build_from_block_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb);

/* ------------------------------------------------------------
 * Build block tree from H^2-matrix
 * ------------------------------------------------------------ */

/** @brief Build an @ref block tree from an @ref h2matrix.
 *
 *  @param G Input @ref h2matrix.
 *  @returns New @ref block tree corresponding to the structure of <tt>G</tt>.
 */
HEADER_PREFIX pblock
build_from_h2matrix_block(pch2matrix G);

/* ------------------------------------------------------------
 * Enumeration by block number
 * ------------------------------------------------------------ */

/** @brief Enumerate @ref h2matrix according to block tree.
 *
 *  The @ref h2matrix submatrices are enumerated in an array of
 *  size <tt>h2->desc</tt>. The enumeration starts with <tt>0</tt> assigned to
 *  the root and then proceeds column-wise starting with <tt>b->sons[0]</tt>
 *  corresponding to the entries <tt>1</tt> to <tt>h2->sons[0]->desc</tt>
 *  in the array and ending with <tt>h2->sons[rsons*csons-1]</tt> corresponding
 *  to the last <tt>b->sons[rsons*csons-1]->desc</tt> entries.
 *
 *  @param h2 Matrix matching the block structure given by <tt>b</tt>.
 *  @returns Array of size <tt>h2->desc</tt> containing pointers to the
 *         @ref h2matrix objects corresponding to descendants of <tt>h2</tt>. */
HEADER_PREFIX ph2matrix *
enumerate_h2matrix( ph2matrix h2);

/* ------------------------------------------------------------
 * Hierarchical iterators
 * ------------------------------------------------------------ */

/** @brief List of @ref h2matrix objects. */
typedef struct _h2matrixlist h2matrixlist;

/** @brief Pointer to @ref h2matrixlist object. */
typedef h2matrixlist *ph2matrixlist;

/** @brief Pointer to constant @ref h2matrixlist object. */
typedef const h2matrixlist *pch2matrixlist;

/** @brief List of @ref h2matrix objects. */
struct _h2matrixlist {
  /** @brief Matrix. */
  ph2matrix G;

  /** @brief Number of the submatrix. */
  uint mname;

  /** @brief Number of row cluster. */
  uint rname;

  /** @brief Number of column cluster. */
  uint cname;

  /** @brief List entry corresponding to father. */
  pch2matrixlist father;

  /** @brief Next item in list. */
  ph2matrixlist next;
};

/** @brief Callback function for @ref h2matrix iterators.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallelization depth.
 *  @param data Additional data. */
typedef void (*h2matrix_callback_t)(ph2matrix G,
				    uint mname, uint rname, uint cname,
				    uint pardepth, void *data);

/** @brief Callback function for @ref h2matrixlist iterators.
 *
 *  @param t Cluster.
 *  @param tname Number of cluster.
 *  @param pardepth Parallelization depth.
 *  @param hl List of @ref h2matrix objects.
 *  @param data Additional data. */
typedef void (*h2matrixlist_callback_t)(pccluster t,
					uint tname,
					uint pardepth,
					pch2matrixlist hl,
					void *data);

/** @brief Iterate through all submatrices of an @ref h2matrix.
 *
 *  If the iterator works with multiple threads, it guarantees that
 *  threads running in parallel call the <tt>pre</tt> and <tt>post</tt>
 *  functions with different row and column clusters.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix block.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallization depth.
 *  @param pre Function called before accessing sons.
 *  @param post Function called after accessing sons.
 *  @param data Additional data passed to callback functions. */
HEADER_PREFIX void
iterate_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		 uint pardepth,
		 h2matrix_callback_t pre, h2matrix_callback_t post,
		 void *data);

/** @brief Iterate through all submatrices of an @ref h2matrix,
 *  collecting all row blocks in an @ref h2matrixlist object.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix block.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallization depth.
 *  @param pre Function called before accessing sons.
 *  @param post Function called after accessing sons.
 *  @param data Additional data passed to callback functions. */
HEADER_PREFIX void
iterate_rowlist_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
			 uint pardepth,
			 h2matrixlist_callback_t pre,
			 h2matrixlist_callback_t post, void *data);

/** @brief Iterate through all submatrices of an @ref h2matrix,
 *  collecting all column blocks in an @ref h2matrixlist object.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix block.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallization depth.
 *  @param pre Function called before accessing sons.
 *  @param post Function called after accessing sons.
 *  @param data Additional data passed to callback functions. */
HEADER_PREFIX void
iterate_collist_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
			 uint pardepth,
			 h2matrixlist_callback_t pre,
			 h2matrixlist_callback_t post, void *data);

/** @brief Iterate through all submatrices of an @ref h2matrix.
 *
 *  If the iterator works with multiple threads, it guarantees that
 *  threads running in parallel call the <tt>pre</tt> and <tt>post</tt>
 *  functions with different row clusters.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix block.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallization depth.
 *  @param pre Function called before accessing sons.
 *  @param post Function called after accessing sons.
 *  @param data Additional data passed to callback functions. */
HEADER_PREFIX void
iterate_byrow_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		       uint pardepth,
		       h2matrix_callback_t pre, h2matrix_callback_t post,
		       void *data);

/** @brief Iterate through all submatrices of an @ref h2matrix.
 *
 *  If the iterator works with multiple threads, it guarantees that
 *  threads running in parallel call the <tt>pre</tt> and <tt>post</tt>
 *  functions with different column clusters.
 *
 *  @param G Matrix.
 *  @param mname Number of matrix block.
 *  @param rname Number of row cluster.
 *  @param cname Number of column cluster.
 *  @param pardepth Parallization depth.
 *  @param pre Function called before accessing sons.
 *  @param post Function called after accessing sons.
 *  @param data Additional data passed to callback functions. */
HEADER_PREFIX void
iterate_bycol_h2matrix(ph2matrix G, uint mname, uint rname, uint cname,
		       uint pardepth,
		       h2matrix_callback_t pre, h2matrix_callback_t post,
		       void *data);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2trans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_h2matrix_avector(field alpha, bool h2trans, pch2matrix h2, pcavector x,
    pavector y);

/** @brief Interaction phase of the matrix-vector multiplication.
 *
 *  Nearfield blocks are added directly
 *  @f$y|_{\hat t} \gets y|_{\hat t} + A|_{\hat t\times\hat s} x|_{\hat s}@f$,
 *  farfield block contributions are accumulated
 *  @f$\hat y_t \gets \hat y_t + S_{t,s} \hat x_s@f$.
 *
 *  Both <tt>xt</tt> and <tt>yt</tt> should be coefficient vectors provided by
 *  @ref new_coeffs_clusterbasis_avector, <tt>xt</tt> is usually initialized by
 *  @ref forward_clusterbasis_avector, while <tt>yt</tt> is typically added to the
 *  result using @ref backward_clusterbasis_avector.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param xt Coefficients @f$(\hat x_s)_{s\in\mathcal{T}_{\mathcal J}}@f$
 *            of the source vector with respect to the
 *            column basis <tt>h2->cb</tt>.
 *  @param yt Coefficients @f$(\hat y_t)_{t\in\mathcal{T}_{\mathcal I}}@f$
 *            of the target vector with respect to the
 *            row basis <tt>h2->rb</tt>. */
HEADER_PREFIX void
fastaddeval_h2matrix_avector(field alpha, pch2matrix h2, pavector xt,
    pavector yt);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_h2matrix_avector(field alpha, pch2matrix h2, pcavector x, pavector y);

/** @brief Interaction phase of the adjoint matrix-vector multiplication.
 *
 *  Nearfield blocks are added directly
 *  @f$y|_{\hat s} \gets y|_{\hat s} + A|_{\hat t\times\hat s}^* x|_{\hat t}@f$,
 *  farfield block contributions are accumulated
 *  @f$\hat y_s \gets \hat y_s + S_{t,s}^* \hat x_t@f$.
 *
 *  Both <tt>xt</tt> and <tt>yt</tt> should be coefficient vectors provided by
 *  @ref new_coeffs_clusterbasis_avector, <tt>xt</tt> is usually initialized by
 *  @ref forward_clusterbasis_avector, while <tt>yt</tt> is typically added to the
 *  result using @ref backward_clusterbasis_avector.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param xt Coefficients @f$(\hat x_t)_{t\in\mathcal{T}_{\mathcal I}}@f$
 *            of the source vector with respect to the
 *            row basis <tt>h2->rb</tt>.
 *  @param yt Coefficients @f$(\hat y_s)_{s\in\mathcal{T}_{\mathcal J}}@f$
 *            of the target vector with respect to the
 *            column basis <tt>h2->cb</tt>. */
HEADER_PREFIX void
fastaddevaltrans_h2matrix_avector(field alpha, pch2matrix h2, pavector xt,
    pavector yt);

/** @brief Adjoint matrix-vector multiplication
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint matrix is multiplied by the source vector @f$x@f$, the
 *  result is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_h2matrix_avector(field alpha, pch2matrix h2, pcavector x,
    pavector y);

/** @brief Symmetric matrix-vector multiplication,
 *  @f$y \gets y + \alpha A x@f$, where @f$A@f$ is assumed to be
 *  self-adjoint and only its lower triangular part is used.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevalsymm_h2matrix_avector(field alpha, pch2matrix h2, pcavector x,
    pavector y);

/* ------------------------------------------------------------
 * Addmul H2-Matrices and Amatrix
 * ------------------------------------------------------------ */

/** @brief Interaction phase of the @ref h2matrix - @ref amatrix multiplication.
 *
 *  If <tt>atrans==false,</tt> nearfield blocks are added directly
 *  @f$C|_{\hat t\times m} \gets C|_{\hat t\times m} + \alpha A|_{\hat t\times\hat s} B|_{\hat s\times m}@f$,
 *  farfield block contributions are accumulated
 *  @f$\widehat Y_t \gets \widehat Y_t + \alpha S_{t,s} \widehat B_s@f$.
 *  Both <tt>bt</tt> and <tt>ct</tt> should be coefficient matrices with the same
 *  number of columns and <tt>cb->ktree</tt> and <tt>rb->ktree</tt> rows, respectively.
 *
 *  If <tt>atrans==true,</tt> nearfield blocks are treated by
 *  @f$C|_{\hat s\times m} \gets C|_{\hat s\times m} + \alpha A|_{\hat t\times\hat s}^* B|_{\hat t\times m}@f$,
 *  farfield block contributions are accumulated
 *  @f$\widehat Y_s \gets \widehat Y_s + \alpha S_{t,s}^* \widehat B_t@f$.
 *  Both <tt>bt</tt> and <tt>ct</tt> should be coefficient matrices with the same
 *  number of columns and <tt>rb->ktree</tt> and <tt>cb->ktree</tt> rows, respectively.
 *
 *  <tt>bt</tt> is usually initialized by @ref forward_clusterbasis_amatrix,
 *  while <tt>ct</tt> is typically added to the result using
 *  @ref backward_clusterbasis_amatrix.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A First source matrix @f$A@f$.
 *  @param Bt Second source matrix, expressed through coefficient
 *            matrices @f$(\widehat X_s)_{s\in\mathcal{T}_{\mathcal J}}@f$
 *            of the second source matrix with respect to the
 *            column basis <tt>h2->cb</tt> if <tt>atrans==false</tt> and through
 *            @f$(\widehat X_t)_{t\in\mathcal{T}_{\mathcal I}}@f$
 *            with respect to the row basis <tt>h2->rb</tt> otherwise.
 *  @param Ct Target matrix, expressed through coefficient matrices
 *            @f$(\widehat Y_t)_{t\in\mathcal{T}_{\mathcal I}}@f$ with
 *            respect to the row basis <tt>h2->rb</tt> if <tt>atrans==false</tt>
 *            and through
 *            @f$(\widehat Y_s)_{s\in\mathcal{T}_{\mathcal J}}@f$ with
 *            respect to the column basis <tt>h2->cb</tt> otherwise. */
HEADER_PREFIX void
fastaddmul_h2matrix_amatrix_amatrix(field alpha, bool atrans, pch2matrix A,
    pcamatrix Bt, pamatrix Ct);

/** @brief Matrix multiplication @f$ C \gets C + \alpha A B @f$,
 *  @f$ C \gets C + \alpha A^* B @f$, @f$ C \gets C + \alpha A B^* @f$ or
 *  @f$ C \gets C + \alpha A^* B^* @f$.
 *
 *  Multiplies an @ref h2matrix by an @ref amatrix.
 *  For the truncated multiplication of two @ref h2matrix objects,
 *  see @ref addmul_h2matrix.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A First source @ref h2matrix @f$A@f$.
 *  @param btrans Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param B Second source @ref amatrix @f$B@f$.
 *  @param C Target @ref amatrix @f$C@f$. */
HEADER_PREFIX void
addmul_h2matrix_amatrix_amatrix(field alpha, bool atrans, pch2matrix A,
    bool btrans, pcamatrix B, pamatrix C);

/** @brief Matrix multiplication @f$ C \gets C + \alpha A B @f$,
 *  @f$ C \gets C + \alpha A^* B @f$, @f$ C \gets C + \alpha A B^* @f$ or
 *  @f$ C \gets C + \alpha A^* B^* @f$.
 *
 *  Multiplies an @ref amatrix by an @ref h2matrix.
 *  For the truncated multiplication of two @ref h2matrix objects,
 *  see @ref addmul_h2matrix.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A First source @ref amatrix @f$A@f$.
 *  @param btrans Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param B Second source @ref h2matrix @f$B@f$.
 *  @param C Target @ref amatrix @f$C@f$. */
HEADER_PREFIX void
addmul_amatrix_h2matrix_amatrix(field alpha, bool atrans, pcamatrix A,
    bool btrans, pch2matrix B, pamatrix C);

/* ------------------------------------------------------------
 * Orthogonal projection
 * ------------------------------------------------------------ */

/** @brief Compute @f$S \gets V_t^* A W_s@f$.
 *
 *  If row and column basis are orthogonal, the result is the coupling
 *  matrix for the optimal approximation of @f$A@f$ in these bases
 *  with respect to the spectral and Frobenius norms.
 *
 *  @param a Source matrix @f$A@f$.
 *  @param rb Row basis @f$(V_t)_{t\in{\mathcal T}_{\mathcal I}}@f$.
 *  @param cb Column basis @f$(W_s)_{s\in{\mathcal T}_{\mathcal J}}@f$.
 *  @param s Target matrix, will be overwritten by
 *         @f$V_t^* A W_s@f$. */
HEADER_PREFIX void
collectdense_h2matrix(pcamatrix a, pcclusterbasis rb, pcclusterbasis cb,
    pamatrix s);

/** @brief Compute the best approximation of a given matrix @f$A@f$
 *  with respect to an @f$\mathcal{H}^2@f$-matrix space.
 *
 *  All nearfield blocks will be copied (and permuted according to
 *  the cluster numbering), all farfield blocks will be approximated
 *  by coupling matrices @f$S_{t,s} = V_t^* A|_{\hat t\times\hat s} W_s@f$.
 *
 *  @remark In order to obtain the best approximation, both the
 *  row basis <tt>h2->rb</tt> and the column basis <tt>h2->cb</tt> should be
 *  orthogonal.
 *
 *  @param h2 Target matrix, will be overwritten by the best
 *         approximation of the source matrix.
 *  @param a Source matrix. */
HEADER_PREFIX void
project_amatrix_h2matrix(ph2matrix h2, pcamatrix a);

/** @brief Compute the best approximation of a given matrix @f$A@f$
 *  with respect to an @f$\mathcal{H}^2@f$-matrix space.
 *
 *  All nearfield blocks will be copied (and permuted according to
 *  the cluster numbering), all farfield blocks
 *  @f$A|_{\hat t\times\hat s} = X Y^*@f$ will be approximated
 *  by coupling matrices @f$S_{t,s} = (V_t^* X) (W_s^* Y)^*@f$.
 *
 *  @remark In order to obtain the best approximation, both the
 *  row basis <tt>h2->rb</tt> and the column basis <tt>h2->cb</tt> should be
 *  orthogonal.
 *
 *  @param h2 Target matrix, will be overwritten by the best
 *         approximation of the source matrix.
 *  @param h Source matrix. */
HEADER_PREFIX void
project_hmatrix_h2matrix(ph2matrix h2, phmatrix h);

/* ------------------------------------------------------------
 * Spectral norm
 * ------------------------------------------------------------ */

/** @brief Approximate the spectral norm @f$\|H\|_2@f$ of a matrix @f$H@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$H^* H@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param H2 @f$\mathcal H^2@f$ matrix @f$H@f$.
 *  @returns Approximation of @f$\|H\|_2@f$. */
HEADER_PREFIX real
norm2_h2matrix(pch2matrix H2);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b @f$\mathcal H^2@f$ matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_h2matrix(pch2matrix a, pch2matrix b);

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
/** @brief Write @ref h2matrix to NetCDF file.
 *
 *  @param G Matrix.
 *  @param name File name. */
HEADER_PREFIX void
write_cdf_h2matrix(pch2matrix G, const char *name);

/** @brief Write @ref h2matrix to part of a NetCDF file.
 *
 *  @param G Matrix.
 *  @param nc_file File handle.
 *  @param prefix Prefix for variable names. */
HEADER_PREFIX void
write_cdfpart_h2matrix(pch2matrix G, int nc_file, const char *prefix);

/** @brief Read @ref h2matrix from NetCDF file.
 *
 *  @param name File name.
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns @ref h2matrix read from file. */
HEADER_PREFIX ph2matrix
read_cdf_h2matrix(const char *name, pclusterbasis rb, pclusterbasis cb);

/** @brief Read @ref h2matrix from part of a NetCDF file.
 *
 *  @param nc_file File handle.
 *  @param prefix Prefix for variable names.
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns @ref h2matrix read from file. */
HEADER_PREFIX ph2matrix
read_cdfpart_h2matrix(int nc_file, const char *prefix,
		      pclusterbasis rb, pclusterbasis cb);

/** @brief Write @ref h2matrix to NetCDF file, including cluster trees
 *  and cluster bases.
 *
 *  @param G Matrix.
 *  @param name File name. */
HEADER_PREFIX void
write_cdfcomplete_h2matrix(pch2matrix G, const char *name);

/** @brief Read @ref h2matrix from NetCDF file, including cluster trees
 *  and cluster bases.
 *
 *  @param name File name.
 *  @returns @ref h2matrix read from file. */
HEADER_PREFIX ph2matrix
read_cdfcomplete_h2matrix(const char *name);
#endif

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
/**
 * @brief Draw a @ref h2matrix to a cairo surface.
 *
 * @param cr Cairo surface to be drawn to.
 * @param G The @ref h2matrix that should be drawn.
 * @param storage Flag that indicates if the storage requirements for every block
 *   should be depicted into the graphic.
 * @param levels Number of levels of the @ref h2matrix that should be drawn.
 *   If @p level == 0 holds, all levels will be drawn.
 */
HEADER_PREFIX void
draw_cairo_h2matrix(cairo_t *cr, pch2matrix G, bool storage, uint levels);
#endif

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

#ifndef H2MATRIX_COMPLETE
#define H2MATRIX_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_h2matrix(pch2matrix) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_h2matrix(pch2matrix) __attribute__ ((const,unused));
#endif

/** @brief Get the number of rows of an @ref h2matrix @f$G@f$.
 *
 *  @param h2 Matrix @f$G@f$.
 *  @return Number of rows of @f$G@f$. */
INLINE_PREFIX uint
getrows_h2matrix(pch2matrix h2)
{
  return h2->rb->t->size;
}

/** @brief Get the number of columns of an @ref h2matrix @f$G@f$.
 *
 *  @param h2 Matrix @f$G@f$.
 *  @returns Number of columns of @f$G@f$. */
INLINE_PREFIX uint
getcols_h2matrix(pch2matrix h2)
{
  return h2->cb->t->size;
}

#endif

/** @} */
