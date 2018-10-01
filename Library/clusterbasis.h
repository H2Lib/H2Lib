
/* ------------------------------------------------------------
 * This is the file "clusterbasis.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file clusterbasis.h
 *  @author Steffen B&ouml;rm
 */

#ifndef CLUSTERBASIS_H
#define CLUSTERBASIS_H

/** @defgroup clusterbasis clusterbasis
 *  @brief Representation of cluster bases for @f$\mathcal{H}^2@f$-matrices.
 *
 *  The @ref clusterbasis class represents a cluster basis
 *  @f$(V_t)_{t\in{\mathcal T}_{\mathcal{I}}}@f$, typically described
 *  by transfer matrices @f$(E_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$
 *  as @f$V_t = \sum_{t'\in\operatorname{sons}(t)} V_{t'} E_{t'}@f$
 *  for non-leaf clusters @f$t@f$.
 *
 *  @{ */

/** @brief Representation of a cluster basis. */
typedef struct _clusterbasis clusterbasis;

/** @brief Pointer to @ref clusterbasis object. */
typedef clusterbasis *pclusterbasis;

/** @brief Pointer to constant @ref clusterbasis object. */
typedef const clusterbasis *pcclusterbasis;

#include "cluster.h"
#include "clusteroperator.h"
#include "uniform.h"
#include "amatrix.h"

/** @brief Representation of a cluster basis. */
struct _clusterbasis {
  /** @brief Corresponding cluster. */
  pccluster t;

  /** @brief Maximal rank */
  uint k;
  /** @brief Sum of ranks in entire subtree below <tt>t</tt> */
  uint ktree;
  /** @brief Maximal rank sum for all branches below <tt>t</tt> */
  uint kbranch;

  /** @brief Leaf matrix @f$V_t@f$ */
  amatrix V;
  /** @brief Transfer matrix @f$E_t@f$ to father */
  amatrix E;

  /** @brief Number of sons, either <tt>t->sons</tt> or zero */
  uint sons;
  /** @brief Pointers to sons */
  pclusterbasis *son;

  /** @brief Weight matrix used to represent total cluster basis */
  pamatrix Z;

  /** @brief References to this cluster basis */
  uint refs;

  /** @brief List of matrices using this basis as row basis */
  puniform rlist;
  /** @brief List of matrices using this basis as column basis */
  puniform clist;
};

/* Now the "clusterbasis" type is completely defined */
#define CLUSTERBASIS_TYPE_COMPLETE

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Initialize a @ref clusterbasis object.
 *
 *  Sets up the components of the object.
 *  If <tt>t</tt> is not a leaf, the array <tt>son</tt> is allocated,
 *  otherwise it is set to null.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusterbasis object. */
HEADER_PREFIX pclusterbasis
init_clusterbasis(pclusterbasis cb, pccluster t);

/** @brief Initialize a @ref clusterbasis object for a leaf.
 *
 *  Sets up the components of the object.
 *  Sets <tt>son</tt> to null.
 *  If <tt>t->sons>0</tt>, this yields a partial cluster basis.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusterbasis object. */
HEADER_PREFIX pclusterbasis
init_leaf_clusterbasis(pclusterbasis cb, pccluster t);

/** @brief Initialize a @ref clusterbasis object representing
 *  a leaf son cluster for a leaf cluster basis.
 *
 *  Sets up the components of the object, making <tt>V</tt> a submatrix
 *  of <tt>src->V</tt>.
 *  Sets <tt>son</tt> to null.
 *  If <tt>t->sons>0,</tt> this yields a partial cluster basis.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param src Source cluster basis.
 *  @param t Corresponding cluster.
 *  @param off Row offset in <tt>src->V,</tt> should be number of elements
 *         of preceding sons.
 *  @returns Initialized @ref clusterbasis object. */
HEADER_PREFIX pclusterbasis
init_sub_clusterbasis(pclusterbasis cb, pclusterbasis src, pccluster t,
		      uint off);

/** @brief Uninitializes a @ref clusterbasis object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  If this @ref clusterbasis references sons, these sons
 *  are unreferences.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be uninitialized.
 *
 *  @param cb Object to be uninitialized. */
HEADER_PREFIX void
uninit_clusterbasis(pclusterbasis cb);

/** @brief Create a new @ref clusterbasis object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusterbasis.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusterbasis) is deleted.
 *
 *  @param t Corresponding cluster.
 *  @return Returns the newly created @ref clusterbasis object.
 */
HEADER_PREFIX pclusterbasis
new_clusterbasis(pccluster t);

/** @brief Create a new @ref clusterbasis object for a leaf.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusterbasis.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusterbasis) is deleted.
 *
 *  @param t Corresponding cluster.
 *  @return Returns the newly created @ref clusterbasis object.
 */
HEADER_PREFIX pclusterbasis
new_leaf_clusterbasis(pccluster t);

/** @brief Delete a @ref clusterbasis object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref clusterbasis references sons, these sons
 *  are unreferenced.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be deleted.
 *
 *  @param cb Object to be deleted. */
HEADER_PREFIX void
del_clusterbasis(pclusterbasis cb);

/* ------------------------------------------------------------
 * Reference counting
 * ------------------------------------------------------------ */

/** @brief Set a pointer to a @ref clusterbasis object, increase its
 *  reference counter, and decrease reference counter of original
 *  pointer target.
 *
 *  @param ptr Pointer to the @ref pclusterbasis variable that will be changed.
 *  @param cb @ref clusterbasis that will be referenced. */
HEADER_PREFIX void
ref_clusterbasis(pclusterbasis *ptr, pclusterbasis cb);

/** @brief Reduce the reference counter of a @ref clusterbasis object.
 *
 *  If the reference counter reaches zero, the object is deleted.
 *
 *  @remark Use @ref ref_clusterbasis with <tt>cb=NULL</tt> instead, since this
 *  guarantees that the pointer is properly invalidated.
 *
 *  @param cb @ref clusterbasis that will be unreferenced. */
HEADER_PREFIX void
unref_clusterbasis(pclusterbasis cb);

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

/** @brief Updates bookkeeping information, e.g., <tt>cb->ktree,</tt> for
 *  a @ref clusterbasis object after its sons have been altered.
 *
 *  @remark This function will not update the sons recursively.
 *  See @ref update_tree_clusterbasis if you really need to do this.
 *  
 *  @param cb @ref clusterbasis that will be updated. */
HEADER_PREFIX void
update_clusterbasis(pclusterbasis cb);

/** @brief Updates bookkeeping information, e.g., <tt>cb->ktree,</tt> for
 *  a @ref clusterbasis object and all its descendants.
 *
 *  @remark If you only want to update <tt>cb</tt> itself, consider
 *  using @ref update_clusterbasis.
 *
 *  @param cb Root of the @ref clusterbasis tree that will be updated. */
HEADER_PREFIX void
update_tree_clusterbasis(pclusterbasis cb);

/** @brief Change the rank of a cluster basis and resize
 *  <tt>cb->V</tt> and <tt>cb->E</tt>, as well as <tt>cb->son[i]->E</tt>
 *  for all sons.
 *
 *  @remark In typical algorithms, the transfer matrices of the sons
 *  and the leaf matrix will subsequently be set to new values.
 *  In order to keep the basis consistent, it is frequently advisable
 *  to copy <tt>cb->E</tt> before calling this routine and updating
 *  the resized matrix afterwards, e.g., by multiplying with the
 *  appropriate basis transformation.
 *
 *  @param cb Cluster basis that will be changed.
 *  @param k New rank, i.e., number of columns of <tt>V</tt> or <tt>cb->son[i]->E</tt>. */
HEADER_PREFIX void
resize_clusterbasis(pclusterbasis cb, uint k);

/** @brief Change the rank of a cluster basis and resize
 *  <tt>cb->V</tt> and <tt>cb->E</tt>, as well as <tt>cb->son[i]->E</tt>
 *  for all sons.
 *
 *  @remark In typical algorithms, the transfer matrices of the sons
 *  and the leaf matrix will subsequently be set to new values.
 *  In order to keep the basis consistent, it is frequently advisable
 *  to copy <tt>cb->E</tt> before calling this routine and updating
 *  the resized matrix afterwards, e.g., by multiplying with the
 *  appropriate basis transformation.
 *
 *  @param cb Cluster basis that will be changed.
 *  @param k New rank, i.e., number of columns of <tt>V</tt> or <tt>cb->son[i]->E</tt>. */
HEADER_PREFIX void
setrank_clusterbasis(pclusterbasis cb, uint k);

/* ------------------------------------------------------------
 * Build clusterbasis based on cluster
 * ------------------------------------------------------------ */

/** @brief Construct a @ref clusterbasis from a cluster tree.
 *
 *  All ranks will be set to zero.
 *
 *  @param t Root cluster.
 *  @returns New root @ref clusterbasis object following the
 *         structure of the cluster tree. */
HEADER_PREFIX pclusterbasis
build_from_cluster_clusterbasis(pccluster t);

/* ------------------------------------------------------------
 * Clone a cluster basis
 * ------------------------------------------------------------ */

/** @brief Create a copy of a @ref clusterbasis.
 *
 *  All sons are copied as well, including all coefficients.
 *
 *  @param cb Source @ref clusterbasis.
 *  @returns Full copy of <tt>cb</tt>. */
HEADER_PREFIX pclusterbasis
clone_clusterbasis(pcclusterbasis cb);

/** @brief Copy the tree structure of a @ref clusterbasis.
 *
 *  All sons are copied as well.
 *  Coefficients are not copied, all ranks of the new @ref clusterbasis
 *  are set to zero.
 *
 *  @param cb Source @ref clusterbasis.
 *  @returns Structure copy of <tt>cb</tt>. */
HEADER_PREFIX pclusterbasis
clonestructure_clusterbasis(pcclusterbasis cb);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get number of active @ref clusterbasis objects.
 *
 *  @returns Number of active @ref clusterbasis objects. */
HEADER_PREFIX uint
getactives_clusterbasis();

/** @brief Get the size of a given @ref clusterbasis object.
 *
 *  Yields the total size of this object and all its descendants.
 *
 *  @param cb @ref clusterbasis root.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_clusterbasis(pcclusterbasis cb);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Delete all weight matrices of the @ref clusterbasis
 *  and its descendants.
 *
 *  @param cb Target @ref clusterbasis object. */
HEADER_PREFIX void
clear_weight_clusterbasis(pclusterbasis cb);

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

/** @brief Hierarchical iterator for a @ref clusterbasis.
 *
 *  Iterate over the cluster basis and all its descendants.
 *  The <tt>pre</tt> function is called for each element before its
 *  descendants are processed, the <tt>post</tt> function is called
 *  afterwards.
 *
 *  <tt>cbname</tt> will take values between <tt>cbname</tt> and <tt>cbname+cb->t->sons-1,</tt>
 *  where <tt>cbname</tt> is interpreted as the index of the root <tt>cb,</tt>
 *  while the indices of the first son start at <tt>cbname+1,</tt> the indices
 *  of the second at <tt>cbname+1+cb->t->son[0]->desc,</tt> and so on.
 *
 *  @param cb Root of the cluster basis tree.
 *  @param cbname Number of the root.
 *  @param pre Function to be called before the descendants of <tt>cb</tt>
 *         are processed.
 *  @param post Function to be called after the descentants of <tt>cb</tt>
 *         have been processed.
 *  @param data Additional data passed on to <tt>pre</tt> and <tt>post</tt>. */
HEADER_PREFIX void
iterate_clusterbasis(pcclusterbasis cb, uint cbname,
		     void (*pre)(pcclusterbasis cb, uint cbname, void *data),
		     void (*post)(pcclusterbasis cb, uint cbname, void *data),
		     void *data);

/** @brief Parallel hierarchical iterator for a @ref clusterbasis.
 *
 *  Iterate over the cluster basis and all its descendants.
 *  The <tt>pre</tt> function is called for each element before its
 *  descendants are processed, the <tt>post</tt> function is called
 *  afterwards.
 *
 *  <tt>cbname</tt> will take values between <tt>cbname</tt> and <tt>cbname+cb->t->sons-1,</tt>
 *  where <tt>cbname</tt> is interpreted as the index of the root <tt>cb,</tt>
 *  while the indices of the first son start at <tt>cbname+1,</tt> the indices
 *  of the second at <tt>cbname+1+cb->t->son[0]->desc,</tt> and so on.
 *
 *  @param cb Root of the cluster basis tree.
 *  @param cbname Number of the root.
 *  @param pardepth Parallelization depth.
 *  @param pre Function to be called before the descendants of <tt>cb</tt>
 *         are processed.
 *  @param post Function to be called after the descentants of <tt>cb</tt>
 *         have been processed.
 *  @param data Additional data passed on to <tt>pre</tt> and <tt>post</tt>. */
HEADER_PREFIX void
iterate_parallel_clusterbasis(pcclusterbasis cb, uint cbname, uint pardepth,
			      void (*pre)(pcclusterbasis cb, uint cbname, void *data),
			      void (*post)(pcclusterbasis cb, uint cbname, void *data),
			      void *data);

/* ------------------------------------------------------------
 * Enumeration by cluster number
 * ------------------------------------------------------------ */

/** @brief Enumerate @ref clusterbasis according to cluster tree.
 *
 *  The @ref clusterbasis elements are enumerated in an array of
 *  size <tt>t->desc</tt>. The enumeration starts with <tt>0</tt> assigned to
 *  the root and then proceeds with <tt>cb->son[0]</tt> corresponding to
 *  the entries <tt>1</tt> to <tt>t->son[0]->desc</tt> in the array and
 *  ending with <tt>cb->son[sons-1]</tt> corresponding to the last
 *  <tt>t->son[sons-1]->desc</tt> entries.
 *
 *  @param t Cluster tree.
 *  @param cb Cluster basis matching the cluster tree given by <tt>t</tt>.
 *  @returns Array of size <tt>t->desc</tt> containing pointers to the
 *         @ref clusterbasis objects corresponding to <tt>cb</tt> and
 *         its descendants. */
HEADER_PREFIX pclusterbasis *
enumerate_clusterbasis(pccluster t, pclusterbasis cb);

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */

/** @brief Create coefficient vector for cluster basis.
 *
 *  Creates a vector of dimension <tt>cb->ktree</tt> to hold the coefficients
 *  computed by @ref forward_clusterbasis_avector and expected by
 *  @ref backward_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @returns Coefficient vector of size <tt>cb->ktree</tt>. */
HEADER_PREFIX pavector
new_coeffs_clusterbasis_avector(pcclusterbasis cb);

/** @brief Forward transformation.
 *
 *  Compute @f$\hat x_t = V_t^* x@f$ for all elements of the cluster basis.
 *  This function also stores the permuted coefficients corresponding to
 *  the leaves of the cluster basis in the result, preparing all the
 *  necessary information for the multiplication phase realized in
 *  @ref fastaddeval_h2matrix_avector and @ref fastaddevaltrans_h2matrix_avector.
 *
 *  If <tt>cb</tt> is not a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>xt</tt> are filled with @f$\hat x_t = V_t^* x@f$.
 *  The following <tt>cb->son[0]->ktree</tt> entries are filled
 *  with the coefficients for the first son, proceeding until the
 *  last <tt>cb->son[sons-1]->ktree</tt> entries are filled with the
 *  coefficients for the last son.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>xt</tt> are also filled with @f$\hat x_t = V_t^* x@f$.
 *  The following <tt>cb->t->size</tt> rows are filled with the
 *  coefficients of the original vector <tt>x</tt> using the cluster
 *  numbering of <tt>cb->t</tt>.
 *
 *  @remark This function accesses the source vector <tt>x</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>x</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in the
 *  original numbering.
 *  In order to work efficiently with subvectors, consider using
 *  @ref forward_nopermutation_clusterbasis_avector .
 *
 *  @param cb Cluster basis.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_clusterbasis_avector(pcclusterbasis cb, pcavector x, pavector xt);

/** @brief Parallel forward transformation.
 *
 *  Parallel version of @ref forward_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients.
 *  @param pardepth Parallelization depth. */
HEADER_PREFIX void
forward_parallel_clusterbasis_avector(pcclusterbasis cb, pcavector x, pavector xt,
			      uint pardepth);

/** @brief Forward transformation for vectors using cluster numbering.
 *
 *  Version of @ref forward_clusterbasis_avector that takes a vector using
 *  cluster numbering.
 *
 *  @param cb Cluster basis.
 *  @param xp Source vector of dimension <tt>cb->t->size</tt> using
 *         cluster numbering corresponding to <tt>cb->t</tt>.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_nopermutation_clusterbasis_avector(pcclusterbasis cb, pcavector xp,
				   pavector xt);

/** @brief Forward transformation without transfer matrices.
 *
 *  Version of @ref forward_clusterbasis_avector that uses the matrices
 *  <tt>cb->V</tt> for all clusters, not only for leaves.
 *
 *  @remark This function is only intended for uniform
 *  hierarchical matrices, not for @f$\mathcal{H}^2@f$-matrices.
 *
 *  @param cb Cluster basis with <tt>cb->V</tt> initialized for all
 *         descendants.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_notransfer_clusterbasis_avector(pcclusterbasis cb, pcavector x, pavector xt);

/** @brief Backward transformation.
 *
 *  Compute @f$y \gets y + V_t \widehat y_t@f$ for all elements
 *  of the cluster basis.
 *  This function also adds permuted coefficients contained in
 *  <tt>yt</tt> to the result vector.
 *  It can be used to add the result of @ref fastaddeval_h2matrix_avector or
 *  @ref fastaddevaltrans_h2matrix_avector to the target vector.
 *
 *  The contents of <tt>yt</tt> are interpreted as in the function
 *  @ref forward_clusterbasis_avector :
 *  if <tt>cb</tt> is not a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>yt</tt> are coefficients to be multiplied by @f$V_t@f$,
 *  followed by <tt>cb->son[0]->ktree</tt> coefficients for the
 *  first son, and <tt>cb->son[i]->ktree</tt> coefficients for the
 *  <tt>i</tt>-th son, until the last son with <tt>i==sons-1</tt>
 *  has been reached.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>yt</tt> are also coefficients to be multiplied by @f$V_t@f$,
 *  followed by <tt>cb->t->size</tt> coefficients in cluster numbering
 *  to be added to appropriate entries of the target vector.
 *
 *  @remark This function accesses the target vector <tt>y</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>y</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in
 *  the original numbering.
 *  In order to work efficiently with subvectors, consider using
 *  @ref backward_nopermutation_clusterbasis_avector .
 *
 *  @param cb Cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector. */
HEADER_PREFIX void
backward_clusterbasis_avector(pcclusterbasis cb, pavector yt, pavector y);

/** @brief Parallel backward transformation.
 *
 *  Parallel version of @ref backward_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector.
 *  @param pardepth Parallelization depth. */
HEADER_PREFIX void
backward_parallel_clusterbasis_avector(pcclusterbasis cb, pavector yt,
			       pavector y, uint pardepth);

/** @brief Backward transformation for vectors in cluster numbering.
 *
 *  Version of @ref backward_clusterbasis_avector that takes a vector using
 *  cluster numbering.
 *
 *  @param cb Cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param yp Target vector of dimension <tt>cb->t->size</tt> using
 *         cluster numbering corresponding to <tt>cb->t</tt>. */
HEADER_PREFIX void
backward_nopermutation_clusterbasis_avector(pcclusterbasis cb, pavector yt,
				    pavector yp);

/** @brief Backward transformation without transfer matrices.
 *
 *  Version of @ref backward_clusterbasis_avector that uses the matrices
 *  <tt>cb->V</tt> for all clusters, not only for leaves.
 *
 *  @remark This function is only intended for uniform
 *  hierarchical matrices, not for @f$\mathcal{H}^2@f$-matrices.
 *
 *  @param cb Cluster basis with <tt>cb->V</tt> initialized for all
 *         descendants.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector. */
HEADER_PREFIX void
backward_notransfer_clusterbasis_avector(pcclusterbasis cb, pavector yt, pavector y);

/* ------------------------------------------------------------
 * Forward and backward transformation for the root only
 * ------------------------------------------------------------ */

/** @brief Compute @f$\hat x_t = V_t^* x@f$.
 *
 *  Similar to @ref forward_nopermutation_clusterbasis_avector , this function
 *  computes @f$\hat x_t = V_t^* x@f$ using the transfer matrices
 *  @f$E_t@f$ and the leaf matrices @f$V_t@f$.
 *
 *  It differs from @ref forward_nopermutation_clusterbasis_avector in that
 *  @f$\hat x_t@f$ is returned only for the root cluster, not for
 *  its descendants.
 *  This allows the function to reduce auxiliary storage requirements
 *  by re-using storage among different branches of the cluster tree.
 *
 *  @param cb Cluster basis.
 *  @param xp Source vector using cluster numbering corresponding
 *         to <tt>cb->t</tt>.
 *  @param xt Target vector of dimension <tt>cb->kbranch</tt>.
 *         Its first <tt>cb->k</tt> rows will be filled by the
 *         transformed coefficients. */
HEADER_PREFIX void
compress_clusterbasis_avector(pcclusterbasis cb, pcavector xp, pavector xt);

/** @brief Add @f$V_t \hat y_t@f$ to target vector @f$y@f$.
 *
 *  Similar to @ref backward_nopermutation_clusterbasis_avector, this function
 *  computes @f$y \gets y + V_t \hat y_t@f$ using the transfer matrices
 *  @f$E_t@f$ and the leaf matrices @f$V_t@f$.
 *
 *  It differs from @ref forward_nopermutation_clusterbasis_avector in that
 *  only @f$V_t \hat y_t@f$ is added, not the contributions of the
 *  descendants.
 *  This allows the function to reduce auxiliary storage requirements
 *  by re-using storage among different branches of the cluster tree.
 *
 *  @param cb Cluster basis.
 *  @param yt Target vector of dimension <tt>cb->kbranch</tt>.
 *         Its first <tt>cb->k</tt> rows contain transformed coefficients.
 *  @param yp Source vector using cluster numbering corresponding
 *         to <tt>cb->t</tt>. */
HEADER_PREFIX void
expand_clusterbasis_avector(pcclusterbasis cb, pavector yt, pavector yp);

/** @brief Compute @f$\widehat{X}_t = V_t^* X@f$.
 *
 *  Version of @ref compress_clusterbasis_avector that works for matrices
 *  instead of vectors, applying the transformation to all columns
 *  simultaneously.
 *
 *  @param cb Cluster basis.
 *  @param Xp Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows.
 *  @param Xt Target vector with <tt>cb->kbranch</tt> rows.
 *         Its first <tt>cb->k</tt> rows will be filled by the
 *         transformed coefficients. */
HEADER_PREFIX void
compress_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp, pamatrix Xt);

/** @brief Compute @f$\widehat{X}_t = V_t^* X@f$.
 *
 *  Parallel version of @ref compress_clusterbasis_amatrix.
 *
 *  @param cb Cluster basis.
 *  @param Xp Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows.
 *  @param Xt Target vector with <tt>cb->kbranch</tt> rows.
 *         Its first <tt>cb->k</tt> rows will be filled by the
 *         transformed coefficients.
 *  @param pardepth Parallelization depth */
void
compress_parallel_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp,
				    pamatrix Xt, uint pardepth);

/* ------------------------------------------------------------
 * Forward and backward transformation for matrices
 * ------------------------------------------------------------ */

/** @brief Matrix forward transformation.
 *
 *  Compute @f$\widehat{X}_t = V_t^* X@f$ for all elements of the
 *  cluster basis.
 *
 *  Matrix version of @ref forward_nopermutation_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param Xp Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows.
 *  @param Xt Target matrix with <tt>cb->ktree</tt> rows, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_clusterbasis_amatrix(pcclusterbasis cb, pcamatrix Xp, pamatrix Xt);

/** @brief Adjoint matrix forward transformation.
 *
 *  Compute @f$\widehat{X}_t = V_t^* X^*@f$ for all elements of the
 *  cluster basis.
 *
 *  Matrix version of @ref forward_nopermutation_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param Xp Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the columns.
 *  @param Xt Target matrix with <tt>cb->ktree</tt> rows, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_clusterbasis_trans_amatrix(pcclusterbasis cb, pcamatrix Xp, pamatrix Xt);

/** @brief Matrix backward transformation.
 *
 *  Compute @f$X \gets X + V_t \widehat{X}_t@f$ for all elements of the
 *  cluster basis.
 *
 *  Matrix version of @ref backward_nopermutation_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param Yt Source matrix with <tt>cb->ktree</tt> rows, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         The matrix will be overwritten by the function.
 *  @param Yp Target matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows. */
HEADER_PREFIX void
backward_clusterbasis_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Yp);

/** @brief Adjoint matrix backward transformation.
 *
 *  Compute @f$X^* \gets X^* + V_t \widehat{X}_t@f$ for all elements of the
 *  cluster basis.
 *
 *  Matrix version of @ref backward_nopermutation_clusterbasis_avector.
 *
 *  @param cb Cluster basis.
 *  @param Yt Source matrix with <tt>cb->ktree</tt> rows, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         The matrix will be overwritten by the function.
 *  @param Yp Target matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the columns. */
HEADER_PREFIX void
backward_clusterbasis_trans_amatrix(pcclusterbasis cb, pamatrix Yt, pamatrix Yp);

/* ------------------------------------------------------------
 * Simple computations
 * ------------------------------------------------------------ */

/** @brief Compute @f$y \gets y + \alpha V_t \hat y_t@f$.
 *
 *  Similar to @ref expand_clusterbasis_avector , but allocates and frees
 *  the required auxiliary storage automatically.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param cb Cluster basis.
 *  @param yc Source vector of dimension <tt>cb->k</tt>.
 *  @param yp Target vector of dimension <tt>cb->t->size</tt> using
 *         cluster numbering corresponding to <tt>cb->t</tt>. */
HEADER_PREFIX void
addeval_clusterbasis_avector(field alpha, pcclusterbasis cb,
		     pcavector yc, pavector yp);

/** @brief Compute @f$\hat y_t \gets \hat y_t + \alpha V_t^* y@f$.
 *
 *  Similar to @ref compress_clusterbasis_avector , but allocates and frees
 *  the required auxiliary storage automatically.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param cb Cluster basis.
 *  @param xp Source vector of dimension <tt>cb->t->size</tt> using
 *         cluster numbering corresponding to <tt>cb->t</tt>.
 *  @param xc Target vector of dimension <tt>cb->k</tt>. */
HEADER_PREFIX void
addevaltrans_clusterbasis_avector(field alpha, pcclusterbasis cb,
			  pcavector xp, pavector xc);

/* ------------------------------------------------------------
 * Orthogonalization
 * ------------------------------------------------------------ */

/** @brief Orthogonalize a cluster basis.
 *
 *  Replace a given cluster basis @f$(V_t)_{t\in\ct}@f$ by an
 *  orthogonal clusterbasis @f$(Q_t)_{t\in\ct}@f$ such that
 *  @f$V_t = Q_t R_t@f$ holds for suitable matrices @f$(R_t)_{t\in\ct}@f$.
 *
 *  @param cb Original cluster basis @f$(V_t)_{t\in\ct}@f$.
 *         Will be overwritten by the orthogonal basis
 *         @f$(Q_t)_{t\in\ct}@f$.
 *  @param co If not null, will be filled with the matrices
 *         @f$(R_t)_{t\in\ct}@f$ describing the basis change
 *         from the old to the new basis.
 *         The structure of this cluster operator has to match
 *         the structure of the cluster basis, you might consider
 *         using @ref build_from_clusterbasis_clusteroperator */
HEADER_PREFIX void
ortho_clusterbasis(pclusterbasis cb, pclusteroperator co);

/** @brief Check whether a cluster basis is orthogonal.
 *
 *  Returns the maximum of @f$\|V_t^* V_t - I\|_F@f$ (for leaf clusters)
 *  and @f$\|\widehat{V}_t^* \widehat{V}_t - I\|_F@f$ (for non-leaf clusters).
 *  If the result is significantly larger than machine eps, the
 *  cluster basis should not be considered orthogonal.
 *
 *  @param cb Cluster basis.
 *  @return Returns a measure for the orthogonality of the @ref clusterbasis
 *    @p cb.
 */
HEADER_PREFIX real
check_ortho_clusterbasis(pcclusterbasis cb);

/** @brief Compute weight matrices for a cluster basis.
 *
 *  Computes matrices @f$(R_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$
 *  for all clusters such that @f$\|V_t \hat x_t\|_2 = \|R_t \hat x_t\|_2@f$
 *  holds for all @f$\hat x_t@f$.
 *
 *  @param cb Cluster basis @f$(V_t)_{t\in{\mathcal T}_{\mathcal{I}}}@f$.
 *  @param co @ref clusteroperator structure matching the
 *         tree of <tt>cb</tt>, will be overwritten with weight
 *         matrices.
 *         A simple way to set up this @ref clusteroperator it
 *         to use @ref build_from_clusterbasis_clusteroperator.
 *  @returns The @ref clusteroperator <tt>co</tt>. */
HEADER_PREFIX pclusteroperator
weight_clusterbasis_clusteroperator(pcclusterbasis cb, pclusteroperator co);

/**
 * @brief Compute weight matrices for a cluster basis using an iterator over the
 * cluster tree.
 *
 * Computes matrices @f$(R_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$
 * for all clusters such that @f$\|V_t \hat x_t\|_2 = \|R_t \hat x_t\|_2@f$
 * holds for all @f$\hat x_t@f$.
 *
 * @param cb Cluster basis @f$(V_t)_{t\in{\mathcal T}_{\mathcal{I}}}@f$.
 * should be computed.
 *
 * @return Returns an array of the weight matrices
 * @f$(R_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$ corresponding to an enumeration
 * of the cluster tree @f$\mathcal{T}_{\mathcal{I}}@f$.
 */
HEADER_PREFIX pamatrix
weight_enum_clusterbasis_clusteroperator(pcclusterbasis cb);

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
/** @brief Write @ref clusterbasis to NetCDF file.
 *
 *  @param cb Cluster basis.
 *  @param name File name. */
HEADER_PREFIX void
write_cdf_clusterbasis(pcclusterbasis cb, const char *name);

/** @brief Write @ref clusterbasis to part of NetCDF file.
 *
 *  @param cb Cluster basis.
 *  @param nc_file File handle.
 *  @param prefix Prefix for variable names. */
HEADER_PREFIX void
write_cdfpart_clusterbasis(pcclusterbasis cb, int nc_file, const char *prefix);

/** @brief Read @ref clusterbasis from NetCdf file.
 *
 *  @param name File name.
 *  @param t Root @ref cluster for cluster basis.
 *  @returns Cluster basis read from file. */
HEADER_PREFIX pclusterbasis
read_cdf_clusterbasis(const char *name, pccluster t);

/** @brief Read @ref clusterbasis from part of a NetCdf file.
 *
 *  @param nc_file File handle.
 *  @param prefix Prefix for variable names.
 *  @param t Root @ref cluster for cluster basis.
 *  @returns Cluster basis read from file. */
HEADER_PREFIX pclusterbasis
read_cdfpart_clusterbasis(int nc_file, const char *prefix, pccluster t);
#endif

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

/* Avoid problems with incomplete type definitions */
#if defined(CLUSTERBASIS_TYPE_COMPLETE) && !defined(CLUSTERBASIS_COMPLETE)
#define CLUSTERBASIS_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX pamatrix
getE_clusterbasis(pclusterbasis) __attribute__ ((unused));
INLINE_PREFIX pamatrix
getV_clusterbasis(pclusterbasis) __attribute__ ((unused));
#endif

/** @brief Get the transfer matrix @f$E@f$ to the father cluster
 * of a @ref clusterbasis.
 *
 * @param cb Cluster basis.
 * @returns Transfer matrix @f$E@f$ */
INLINE_PREFIX pamatrix
getE_clusterbasis(pclusterbasis cb)
{
  return &cb->E;
}

/** @brief Get the leaf matrix @f$V@f$ of a @ref clusterbasis.
 *
 * @param cb Cluster basis.
 * @returns Leaf matrix @f$V@f$ */
INLINE_PREFIX pamatrix
getV_clusterbasis(pclusterbasis cb)
{
  return &cb->V;
}

#endif

/** @} */

