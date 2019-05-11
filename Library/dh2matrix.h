
/* ------------------------------------------------------------
 * This is the file "dh2matrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file dh2matrix.h
 *  @author Steffen B&ouml;rm */

#ifndef DH2MATRIX_H
#define DH2MATRIX_H

/** @defgroup dh2matrix dh2matrix
 *  @brief @f$\mathcal{DH}^2@f$-matrices.
 *  @{ */

/** @brief @f$\mathcal{DH}^2@f$-matrix. */
typedef struct _dh2matrix dh2matrix;

/** @brief Pointer to a @ref dh2matrix */
typedef dh2matrix *pdh2matrix;

/** @brief Pointer to a constant @ref dh2matrix */
typedef const dh2matrix *pcdh2matrix;

#include "settings.h"
#include "amatrix.h"
#include "krylov.h"
#include "h2matrix.h"
#include "dblock.h"
#include "duniform.h"
#include "dclusterbasis.h"
#include "dclusteroperator.h"

/** @brief Tree structure representing a @f$\mathcal{DH}^2@f$-matrix.
 *
 *  If <tt>u</tt> is not null, we are dealing with an admissible
 *  submatrix represented as @f$V_{tc} S_{ts} W_{sc}^*@f$, where
 *  the row basis @f$V_{tc}@f$ and the column basis @f$W_{sc}@f$ are
 *  given by <tt>rb</tt> and <tt>cb</tt>, respectively, and
 *  @f$S_{ts}@f$ is contained in <tt>u</tt>.
 *
 *  If <tt>f</tt> is not null, we are dealing with an inadmissible
 *  leaf represented directly by <tt>f</tt>.
 *
 *  If <tt>son</tt> is not null, the matrix is subdivided into
 *  <tt>rsons</tt> row sons and <tt>csons</tt> column sons.
 *  The pointer to the submatrix in the <tt>i</tt>-th row and
 *  the <tt>j</tt>-th column can be found <tt>son[i+j*rsons]</tt>. */
struct _dh2matrix {
  /** @brief Row cluster basis */
  pdclusterbasis rb;
  /** @brief Column cluster basis */
  pdclusterbasis cb;

  /** @brief Uniform matrix (if farfield leaf) */
  pduniform u;
  /** @brief Standard matrix (if nearfield leaf) */
  pamatrix f;

  /** @brief Son matrices (if subdivided) */
  pdh2matrix *son;
  /** @brief Block rows */
  uint rsons;
  /** @brief Block columns */
  uint csons;

  /** @brief Number of descendants */
  uint desc;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new empty @ref dh2matrix object.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New object. */
HEADER_PREFIX pdh2matrix
new_dh2matrix(pdclusterbasis rb, pdclusterbasis cb);

/** @brief Create a new @ref dh2matrix object representing an
 *    admissible leaf.
 *
 *  @param rb Row cluster basis.
 *  @param rd Row direction (for directional cluster basis <tt>rb</tt>).
 *  @param cb Column cluster basis.
 *  @param cd Column direction (for directional cluster basis <tt>cb</tt>).
 *  @returns New object. */
HEADER_PREFIX pdh2matrix
new_uniform_dh2matrix(pdclusterbasis rb, uint rd, pdclusterbasis cb, uint cd);

/** @brief Create a new @ref dh2matrix object representing an
 *    inadmissible leaf.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New object. */
HEADER_PREFIX pdh2matrix
new_full_dh2matrix(pdclusterbasis rb, pdclusterbasis cb);

/** @brief Create a new @ref dh2matrix object representing a
 *    subdivided matrix.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @param rsons Number of row sons.
 *  @param csons Number of column sons.
 *  @returns New object. */
HEADER_PREFIX pdh2matrix
new_super_dh2matrix(pdclusterbasis rb, pdclusterbasis cb,
		    uint rsons, uint csons);

/** @brief Update internal data structures after the sons
 *    have been modified, e.g., to recompute <tt>desc</tt>.
 *
 *  @param h2 Matrix. */
HEADER_PREFIX void
update_dh2matrix(pdh2matrix h2);

/** @brief Delete a @ref dh2matrix object.
 *
 *  @param h2 Object to be deleted. */
HEADER_PREFIX void
del_dh2matrix(pdh2matrix h2);

/** @brief Builds a new @ref dh2matrix with the same structure as <tt>dh2</tt> 
 *    and @ref dclusterbasis <tt>rb</tt> as directional row cluster basis and  
 *    <tt>cb</tt> for the columns. 
 * 
 *  All coupling and nearfield matrices from <tt>dh2</tt> will be copied.
 *
 *  @param dh2 @ref dh2matrix which should be cloned.
 *  @param rb Directional row cluster basis for the clone.
 *  @param cb Directional column cluster basis for the clone.
 *  @returns New @ref dh2matrix with the same structure and submatricies as
 *    <tt>dh2</tt>. */
HEADER_PREFIX pdh2matrix
clone_dh2matrix(pcdh2matrix dh2, pdclusterbasis rb, pdclusterbasis cb);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Compute the storage size of a @ref dh2matrix, without
 *    the cluster bases.
 *
 *  @param h2 Matrix.
 *  @returns Storage size in bytes. */
HEADER_PREFIX size_t
getsize_dh2matrix(pcdh2matrix h2);

/** @brief Compute the storage size of the inadmissible leaves of a
 *    @ref dh2matrix.
 *
 *  @param h2 Matrix.
 *  @returns Nearfield storage size in bytes. */
HEADER_PREFIX size_t
getnearsize_dh2matrix(pcdh2matrix h2);

/** @brief Compute the storage size of the admissible leaves of a
 *    @ref dh2matrix.
 *
 *  @param h2 Matrix.
 *  @returns Farfield storage size in bytes. */
HEADER_PREFIX size_t
getfarsize_dh2matrix(pcdh2matrix h2);

/** @brief Compute the storage size of a @ref dh2matrix,
 *    including the cluster bases.
 *
 *  @param h2 Matrix.
 *  @returns Storage size in bytes. */
HEADER_PREFIX size_t
gettotalsize_dh2matrix(pcdh2matrix h2);

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_dh2matrix(pcdh2matrix h2) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_dh2matrix(pcdh2matrix h2) __attribute__ ((const,unused));
#endif

/** @brief Return the number of rows of a @ref dh2matrix.
 *
 *  @param h2 Matrix.
 *  @returns Number of rows. */
INLINE_PREFIX uint
getrows_dh2matrix(pcdh2matrix h2)
{
  return h2->rb->t->size;
}

/** @brief Return the number of columns of a @ref dh2matrix.
 *
 *  @param h2 Matrix
 *  @returns Number of columns */
INLINE_PREFIX uint
getcols_dh2matrix(pcdh2matrix h2)
{
  return h2->cb->t->size;
}

/* ------------------------------------------------------------
 * Build directional H^2-matrix based on directional block tree
 * ------------------------------------------------------------ */

/** @brief Construct a @ref dh2matrix based on a @ref dblock tree
 *    and cluster bases.
 *
 *  @param b Directional block tree.
 *  @param rb Row basis.
 *  @param cb Column basis.
 *  @returns New @ref dh2matrix. */
HEADER_PREFIX pdh2matrix
buildfromblock_dh2matrix(pcdblock b, pdclusterbasis rb, pdclusterbasis cb);

/* ------------------------------------------------------------
 * Copy nearfield matrices from a given dense matrix
 * ------------------------------------------------------------ */

/** @brief Copy nearfield matrices from a given matrix.
 *
 *  If we approximate a given dense matrix @f$G@f$ by a
 *  @f$\mathcal{DH}^2@f$-matrix, we can copy the nearfield
 *  matrices directly from @f$G@f$ instead of recomputing them.
 *
 *  @param G Source matrix @f$G@f$.
 *  @param Gh Target @f$\mathcal{DH}^2@f$-matrix. */
HEADER_PREFIX void
copynear_dh2matrix(pcamatrix G, pdh2matrix Gh);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Fast matrix-vector multiplication.
 *
 *  This function handles the coupling phase of the matrix-vector
 *  multiplication. It is assumed that the transformed input vector
 *  is provided in @f$xt@f$, and the results are added to the
 *  output vector @f$yt@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param xt Transformed input vector @f$\widehat{x}@f$ as computed
 *    by @ref forward_dclusterbasis.
 *  @param yt Transformed output vector @f$\widehat{y}@f$ as required
 *    by @ref backward_dclusterbasis. */
HEADER_PREFIX void
fastaddeval_dh2matrix_avector(field alpha, pcdh2matrix h2,
			      pcavector xt, pavector yt);

/** @brief Fast adjoint matrix-vector multiplication.
 *
 *  This function handles the coupling phase of the matrix-vector
 *  multiplication. It is assumed that the transformed input vector
 *  is provided in @f$xt@f$, and the results are added to the
 *  output vector @f$yt@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param xt Transformed input vector @f$\widehat{x}@f$ as computed
 *    by @ref forward_dclusterbasis.
 *  @param yt Transformed output vector @f$\widehat{y}@f$ as required
 *    by @ref backward_dclusterbasis. */
HEADER_PREFIX void
fastaddevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector xt, pavector yt);

/** @brief Matrix-vector multiplication @f$y \gets y + \alpha G x@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
addeval_dh2matrix_avector(field alpha, pcdh2matrix h2,
			  pcavector x, pavector y);

/** @brief Adjoint matrix-vector multiplication @f$y \gets y + \alpha G^* x@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2,
			       pcavector x, pavector y);

/** @brief Matrix-vector multiplication @f$y \gets y + \alpha G x@f$
 *    or @f$y \gets y + \alpha G^* x@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2trans Set if @f$G^*@f$ is to be used instead of @f$G@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
mvm_dh2matrix_avector(field alpha, bool h2trans, pcdh2matrix h2,
		      pcavector x, pavector y);

/* ------------------------------------------------------------
 * Matrix-vector multiplication based on the
 * dclusterbasis parallel iterator
 * ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication @f$y \gets y + \alpha G x@f$,
 *    parallelized version.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
addeval_parallel_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector x, pavector y);

/** @brief Adjoint matrix-vector multiplication @f$y \gets y + \alpha G^* x@f$,
 *    parallelized version.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_parallel_dh2matrix_avector(field alpha, pcdh2matrix h2,
					pcavector x, pavector y);

/* ------------------------------------------------------------
 * Slow direct matrix-vector multiplication,
 * for debugging purposes
 * ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication @f$y \gets y + \alpha G x@f$,
 *    slow version for debugging, all submatrices are handled
 *    independently.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
slowaddeval_dh2matrix_avector(field alpha, pcdh2matrix h2,
			      pcavector x, pavector y);

/** @brief Adjoint matrix-vector multiplication @f$y \gets y + \alpha G^* x@f$,
 *    slow version for debugging, all submatrices are handled
 *    independently.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 @f$\mathcal{DH}^2@f$-matrix @f$G@f$.
 *  @param x Input vector @f$x@f$.
 *  @param y Output vector @f$y@f$. */
HEADER_PREFIX void
slowaddevaltrans_dh2matrix_avector(field alpha, pcdh2matrix h2,
				   pcavector x, pavector y);


/* ------------------------------------------------------------
 * Conversion to a full matrix
 * ------------------------------------------------------------ */

/** @brief Add a @f$\mathcal{DH}^2@f$-matrix to a standard matrix,
 *    @f$G \gets G + \alpha \widetilde{G}@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$\widetilde{G}@f$.
 *  @param G Target matrix @f$G@f$. */
HEADER_PREFIX void
expand_dh2matrix(field alpha, pcdh2matrix h2, pamatrix G);

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

/** @brief Enumerate all submatrices of a @ref dh2matrix tree.
 *
 *  The root is placed at the beginning of the array, followed
 *  by the enumerations of the sons in column-major order.
 *
 *  @param h2 Matrix.
 *  @returns Array of <tt>h2->desc</tt> pointers to all submatrices
 *    of <tt>h2</tt>. */
HEADER_PREFIX pdh2matrix *
enumerate_dh2matrix(pdh2matrix h2);

/** @brief Apply functions to all submatrices of a @ref dh2matrix tree.
 *
 *  The function <tt>pre</tt> is called before the sons are processed,
 *  the function <tt>post</tt> is called afterwards.
 *
 *  @param G Matrix.
 *  @param mname Number of the submatrix, matching the order used
 *    by @ref enumerate_dh2matrix.
 *  @param rname Number of the row cluster, matching the order used
 *    by @ref enumerate_cluster and @ref enumerate_dclusterbasis.
 *  @param cname Number of the column cluster, matching the order used
 *    by @ref enumerate_cluster and @ref enumerate_dclusterbasis.
 *  @param pardepth Maximal depth for recursive parallelization.
 *  @param pre Function to be called before visiting the sons.
 *  @param post Function to be called after the sons have been visited.
 *  @param data Additional data passed on to <tt>pre</tt> and <tt>post</tt>.*/
HEADER_PREFIX void
iterate_dh2matrix(pdh2matrix G, uint mname, uint rname, uint cname,
		  uint pardepth,
		  void (*pre)(pdh2matrix G, uint mname,
			      uint rname, uint cname, uint pardepth,
			      void *data),
		  void (*post)(pdh2matrix G, uint mname,
			       uint rname, uint cname, uint pardepth,
			       void *data), void *data);

/* ------------------------------------------------------------
 * Matrix norms
 * ------------------------------------------------------------ */

/** @brief Approximate the spectral norm @f$\|H\|_2@f$ of a matrix @f$H@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$H^* H@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param DH2 Directional @f$\mathcal H^2@f$ matrix @f$H@f$.
 *  @returns Approximation of @f$\|H\|_2@f$. */
HEADER_PREFIX real
norm2_dh2matrix(pcdh2matrix DH2);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a D@f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b D@f$\mathcal H^2@f$ matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_dh2matrix(pcdh2matrix a, pcdh2matrix b);

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
/** @brief Draw a @ref dh2matrix tree to a cairo surface.
 *
 *  @param cr Cairo surface to be drawn to.
 *  @param G Matrix to be drawn.
 *  @param storage Set if the storage requirements of the leaf
 *    submatrices are to be displayed.
 *  @param ranks Flag that indicates if the size of the row and column cluster 
 *    should be depicted into the graphic.   
 *  @param levels Number of levels of the @ref dh2matrix that should be drawn.
 *    For <tt>level==0</tt>, all levels will be drawn. */
HEADER_PREFIX void
draw_cairo_dh2matrix(cairo_t * cr, pcdh2matrix G, bool storage,
		     bool ranks, uint levels);
#endif

/* -----------------------------------------------------------
 * Resize, projection and recompression
 * ----------------------------------------------------------- */

/** 
 *  @brief Computes new coupling matrices @f$S_b@f$ after
 *         the directional cluster basis has been orthogonalized.
 *         
 *  Computes the new coupling matrices @f$S_b@f$ for the orthogonal
 *  directional cluster basis for all admissible directional blocks 
 *  @f$ b = (t,s) @f$ with direction @f$ \iota_{t} @f$
 *  through @f$ R_{t,\iota_{t}} S_b R^{*}_{s,\iota_{t}} @f$, 
 *  where @f$R_{t, \iota_{t}} @f$ is the matrix describing the 
 *  basis change of @f$ V_{t, \iota_{t}} @f$. 
 *  This basis change matrices are computed during the orthogonalization and
 *  saved in a @ref dclusteroperator object.
 *  
 *  @remark Should be called after orthogonalize a given directional 
 *         cluster basis with @ref ortho_dclusterbasis.
 *  @param A Directional @f$\mathcal{H}^2@f$-matrix including the coupling matrices.
 *  @param ro Directional cluster operator for the directional row cluster basis.
 *  @param co Directional cluster operator for the directional column cluster basis.       
 */

HEADER_PREFIX void
resize_coupling_dh2matrix(pdh2matrix A, pdclusteroperator ro, pdclusteroperator co);

/** 
 *  @brief Builds a new @f$\mathcal{DH}^2@f$-matrix approximation
 *         of a given @f$\mathcal{DH}^2@f$-matrix with new cluster bases.
 *         
 *  The structure and nearfield matrices of the old @f$\mathcal{DH}^2@f$-matrix
 *  are copied, new coupling matrices are computed with the basis change
 *  matrices saved in the directional cluster operator objects.
 *
 *  @param dh2 Original @f$\mathcal{DH}^2@f$-matrix.
 *  @param rb New directional row basis.
 *  @param ro @ref dclusteroperator describing the basis change from 
 *         old row basis <tt>dh2->rb</tt> to the new basis <tt>rb</tt>.
 *  @param cb New directional column basis.
 *  @param co @ref dclusteroperator describing the basis change from 
 *         old column basis <tt>dh2->cb</tt> to the new basis <tt>cb</tt>.
 *  @returns The new @f$\mathcal{DH}^2@f$-matrix approximation. 
 */

HEADER_PREFIX pdh2matrix
build_projected_dh2matrix(pcdh2matrix dh2, pdclusterbasis rb, pdclusterbasis cb, pdclusteroperator ro, pdclusteroperator co);

/** 
 *  @brief Computes a @f$\mathcal{DH}^2@f$-matrix approximation of 
 *  a @f$\mathcal{DH}^2@f$-matrix.
 *
 *  @param G Source @ref dh2matrix @f$G@f$.
 *  @param rbortho Set if the original row basis is orthogonal,
 *         this allows the algorithm to avoid computing row weights.
 *  @param cbortho Set if the original column basis is orthogonal,
 *         this allows the algorithm to avoid computing column weights.
 *  @param tm Truncation mode.
 *  @param eps Truncation accuracy.
 *  @returns @f$\mathcal{DH}^2@f$-matrix approximation of @f$G@f$. 
 */


HEADER_PREFIX pdh2matrix
compress_dh2matrix_dh2matrix(pcdh2matrix G, bool rbortho, bool cbortho, pctruncmode tm, real eps);

/**
 *  @brief Computes truncated directional cluster basis and save basis change matrices.
 * 
 *  For every cluster @f$ t @f$ and direction @f$ \iota @f$ a singular value 
 *  decomposition of the directional cluster basis multiplied by a weight matrix 
 *  stored in <tt>co</tt> will be made. 
 *  Then the smallest possible rank @f$ k_new @f$ for the desired accuracy and 
 *  truncation mode is selected.
 *  The new directional cluster basis and transfer matrices for this rank @f$ k_new @f$ 
 *  are computed and saved in the @ref dclusterbasis object <tt>cnew</tt>. 
 *  Matrices describing the basis change are stored in <tt>bco</tt>.
 * 
 *  @param cold Old @ref dclusterbasis object.
 *  @param cnew @ref dclusterbasis object for the truncated cluster basis.
 *  @param co @ref dclusteroperator object including weight matrices.
 *  @param bco Directional cluster operator for saving basis change matrices.
 *  @param tm @ref truncmode object for applied mode of truncation.
 *  @param eps Real value used as accuracy for the truncation.
 */

HEADER_PREFIX void
truncate_dclusterbasis(pdclusterbasis cold, pdclusterbasis cnew, pdclusteroperator co, pdclusteroperator bco, pctruncmode tm, real eps);

/** @} */

#endif
