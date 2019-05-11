
/* ------------------------------------------------------------
 * This is the file "dclusterbasis.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file dclusterbasis.h
 *  @author Steffen B&ouml;rm
 */

#ifndef DCLUSTERBASIS_H
#define DCLUSTERBASIS_H

/** @defgroup dclusterbasis dclusterbasis
 *  @brief Directional cluster bases for directional @f$\mathcal{H}^2@f$-matrices.
 *  @{ */

/** @brief Representation of a directional cluster basis */
typedef struct _dclusterbasis dclusterbasis;

/** @brief Pointer to a @ref dclusterbasis object */
typedef dclusterbasis *pdclusterbasis;

/** @brief Pointer to a constant @ref dclusterbasis object */
typedef const dclusterbasis *pcdclusterbasis;

#include "dcluster.h"
#include "dblock.h"
#include "amatrix.h"
#include "dclusteroperator.h"

/** @brief Representation of a directional cluster basis.
 *
 *  A directional cluster basis is a family
 *  @f$(V_{t,\iota})_{t\in\ctI,\iota\in D_t}@f$
 *  of matrices that assigns each cluster @f$t\in\ctI@f$ and each direction
 *  @f$\iota\in D_t@f$ corresponding to this cluster a basis of a subspace
 *  of @f$\bbbr^{\hat t}@f$ spanned by the columns of @f$V_{t,\iota}@f$.
 *
 *  For the sake of efficiency, cluster bases are usually stored
 *  in nested representation, i.e., if @f$t@f$ has sons, @f$V_{t,\iota}@f$
 *  is represented @e implicitly by @f$V_{t,\iota}|_{\hat t'}=V_{t',\kappa} E_{t',\iota}@f$
 *  using transfer matrices @f$E_{t',\iota}@f$, where @f$\kappa \in D_{t'}@f$ is
 *  a suitable direction of the son @f$t'@f$. */
struct _dclusterbasis
{
  /** @brief Corresponding directional cluster */
  pcdcluster t;

  /** @brief Number of directions, matches <tt>t->directions</tt> if
   *  <tt>t->direction > 0</tt>, otherwise equals one. */
  uint directions;

  /** @brief Ranks, i.e., number of columns of @f$V_{t,\iota}@f$,
   *    type <tt>uint k[directions]</tt> */
  uint *k;

  /** @brief Partial sums of direction ranks for @ref forward_dclusterbasis
   *    and @ref backward_dclusterbasis,
   *    <tt>koff[iota]=</tt>@f$\sum_{\kappa=0}^{\iota-1} k_\kappa@f$,
   *    type <tt>uint koff[directions+1]</tt> */
  uint *koff;

  /** @brief Sum of ranks in the entire subtree */
  uint ktree;

  /** @brief Maximum of rank sums along all branches in the subtree */
  uint kbranch;
        
  /** @brief Matrices @f$V_{t,\iota}@f$, only stored if <tt>t->sons==0</tt>,
   *    type <tt>amatrix V[directions]</tt> */
  pamatrix V;

  /** @brief Transfer matrices @f$E_{t',\iota}@f$ for all sons,
   *    type <tt>amatrix E[sons][directions]</tt> */
  pamatrix *E;

  /** @brief Number of sons */
  uint sons;
  /** @brief Pointers to sons */
  pdclusterbasis *son;
  /** @brief Son directions corresponding to this cluster's directions,
   *    type <tt>uint dirson[sons][directions]</tt> */
  uint **dirson;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Initialize a @ref dclusterbasis object.
 *
 *  Sets up all components of the dclusterbasis object.
 *  If <tt>t</tt> is not a leaf, the array <tt>son</tt> is allocated,
 *  otherwise it is set to null.
 *
 *  @attention To initialize the matrices for every @f$t@f$ and direction
 *  @f$\iota\in D_t@f$ the ranks have to be find.
 *  Therefore the function @ref findranks_dclusterbasis and
 *  @ref initmatrices_dclusterbasis exists.
 *  
 *  @remark Should always be matched by a call to @ref uninit_dclusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref dclusterbasis object. */

HEADER_PREFIX pdclusterbasis
init_dclusterbasis(pdclusterbasis cb, pcdcluster t);

/** @brief Uninitializes a @ref dclusterbasis object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param cb Object to be uninitialized. */

HEADER_PREFIX void
uninit_dclusterbasis(pdclusterbasis cb);

/** @brief Create a new @ref dclusterbasis object.
 *
 *  Allocates storage for a new directional cluster basis
 *  and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_dclusterbasis.
 *
 *  @param t Corresponding directional cluster.
 *  @return Returns the newly created @ref dclusterbasis object. */

HEADER_PREFIX pdclusterbasis
new_dclusterbasis(pcdcluster t);

/** @brief Delete a @ref dclusterbasis object.
 *
 *  Releases the storage corresponding to the directional 
 *  cluster basis.
 *
 *  @param cb Object to be deleted. */

HEADER_PREFIX void
del_dclusterbasis(pdclusterbasis cb);

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

/** @brief Updates bookkeeping information, e.g., <tt>cb->ktree,</tt> for
 *  a @ref dclusterbasis object after its sons have been altered.
 *
 *  @remark This function will not update the sons recursively!
 *  
 *  @param cb @ref dclusterbasis that will be updated. */

HEADER_PREFIX void
update_dclusterbasis(pdclusterbasis cb);

/** @brief Change the rank of a directional cluster basis and resize
 *  <tt>cb->V[iota]</tt>, while for all sons <tt>i</tt> the transfer
 *  matrices <tt>cb->E[i][iota]</tt> is resized.
 *
 *  @param cb Directional cluster basis that will be changed.
 *  @param iota Direction whose rank is changed.
 *  @param k New rank, i.e., number of columns of <tt>V</tt> or <tt>cb->E[i][iota]</tt>. */

HEADER_PREFIX void
setrank_dclusterbasis(pdclusterbasis cb, uint iota, uint k);

/** @brief Initialize the matrices for a @ref dclusterbasis object.
 * 
 *  Resize the matrices for the whole directional cluster basis
 *  so that all cluster basis @f$V@f$ and transfer matrices
 *  @f$E@f$ are initialized with the desired size and are ready
 *  to be filled.
 *    
 *  @param cb Directional cluster basis, which matrices will be
 *  initialized. */

HEADER_PREFIX void
initmatrices_dclusterbasis(pdclusterbasis cb);

/** @brief Determine the ranks of the directional cluster basis.
 *
 *  @attention This function does not initialize the transfer
 *  and leaf matrices. Functions like @ref initmatrices_dclusterbasis
 *  take care of this part of the initialization. 
 *  @param rank Values for the rank.
 *  @param b Corresponding @ref dblock object.
 *  @param rb Directional row cluster basis.
 *  @param cb Directional column cluster basis. */

HEADER_PREFIX void 
  findranks_dclusterbasis(uint rank, pcdblock b, pdclusterbasis rb,
                          pdclusterbasis cb); 


/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */
/** @brief Get number of active @ref dclusterbasis objects.
 *
 *  @returns Number of active @ref dclusterbasis objects. */

HEADER_PREFIX uint
getactives_dclusterbasis();

/** @brief Get the size of a given @ref dclusterbasis object.
 *
 *  Yields the total size of this object and all its descendants.
 *
 *  @param cb @ref dclusterbasis root.
 *  @returns Size of allocated storage in bytes. */

HEADER_PREFIX size_t
getsize_dclusterbasis(pcdclusterbasis cb);

/** @brief Get the size of a given @ref dclusterbasis object, without 
 * considering its sons.
 *
 *  Yields the total size of a directional cluster basis without 
 *  taking all its descendants into account.
 *
 *  @param cb @ref dclusterbasis root.
 *  @returns Size of allocated storage in bytes. */

HEADER_PREFIX size_t
getsize_nonrecursive_dclusterbasis(pcdclusterbasis cb);

/** @brief Get the maximum rank of a given @ref dclusterbasis object.
 *
 *  Yields the maximum rank for all directions and descendants
 *  of a directional cluster basis.
 *
 *  @param cb @ref dclusterbasis root.
 *  @returns Maximum rank. */

HEADER_PREFIX uint
getmaxrank_dclusterbasis(pcdclusterbasis cb);

/** @brief Count the number of all used directions for
 *  a given @ref dclusterbasis object.
 *
 *  @param cb @ref dclusterbasis root.
 *  @returns Number of all used directions. */

HEADER_PREFIX uint
getactivedirections_dclusterbasis(pcdclusterbasis cb);

/** @brief Print the tree structure of a given 
 *  @ref dclusterbasis object.
 *  
 *  Prints the structure of the given directional
 *  cluster basis and for every basis and directions
 *  his rank, <tt>cb->ktree</tt> and the size of the
 *  corresponding cluster. 
 *  
 *  @param cb Directional cluster basis to be printed.*/

HEADER_PREFIX void
print_tree_dclusterbasis(pcdclusterbasis cb);

/* ------------------------------------------------------------
 * Build cluster basis based on directional cluster tree
 * ------------------------------------------------------------ */

/** @brief Construct a @ref dclusterbasis from a 
 *  directional cluster tree.
 *  
 *  Only the structure is for the @ref dclusterbasis object
 *  is build, no ranks are determined or matrices initialized.
 *
 *  @param t Root directional cluster.
 *  @returns New root @ref dclusterbasis object following the
 *         structure of the directional cluster tree. */

HEADER_PREFIX pdclusterbasis
buildfromdcluster_dclusterbasis(pcdcluster t);

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */
/** @brief Create coefficient vector for a directional cluster basis.
 *
 *  Creates a vector of dimension <tt>cb->ktree</tt> to hold the coefficients
 *  computed by @ref forward_dclusterbasis and expected by
 *  @ref backward_dclusterbasis.
 *
 *  @param cb Directional custer basis.
 *  @returns Coefficient vector of size <tt>cb->ktree</tt>. */

HEADER_PREFIX pavector
newcoeffs_dclusterbasis(pcdclusterbasis cb);

/** @brief Forward transformation.
 *
 *  Compute @f$\hat x_{t \iota} = V_{t \iota}^* x@f$ for all elements of 
 *  the cluster basis.
 *  This function also stores the permuted coefficients corresponding to
 *  the leaves of the cluster basis in the result, preparing all the
 *  necessary information for the multiplication phase realized in
 *  @ref fastaddeval_dh2matrix_avector and @ref fastaddevaltrans_dh2matrix_avector.
 *
 *  If <tt>cb</tt> is not a leaf, the first <tt>cb->koff[cb->directions]</tt> rows of
 *  <tt>xt</tt> are successively filled with @f$\hat x_{t \iota} = V_{t \iota}^* x@f$.
 *  The following <tt>cb->son[0]->ktree</tt> entries are filled
 *  with the coefficients for the first son, proceeding until the
 *  last <tt>cb->son[sons-1]->ktree</tt> entries are filled with the
 *  coefficients for the last son.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->koff[cb->directions]</tt> rows of
 *  <tt>xt</tt> are also successively filled with @f$\hat x_{t \iota} = V_{t \iota}^* x@f$.
 *  The following <tt>cb->t->size</tt> rows are filled with the
 *  coefficients of the original vector <tt>x</tt> using the cluster
 *  numbering of <tt>cb->t</tt>.
 *
 *  @remark This function accesses the source vector <tt>x</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>x</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in the
 *  original numbering..
 *
 *  @param cb Directional cluster basis.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */

HEADER_PREFIX void
forward_dclusterbasis(pcdclusterbasis cb, pcavector x, pavector xt);

/** @brief Backward transformation.
 *
 *  Compute @f$y \gets y + V_{t \iota} \widehat y_{t \iota}@f$ for all elements
 *  of the cluster basis.
 *  This function also adds permuted coefficients contained in
 *  <tt>yt</tt> to the result vector.
 *  It can be used to add the result of @ref fastaddeval_dh2matrix_avector or
 *  @ref fastaddevaltrans_dh2matrix_avector to the target vector.
 *
 *  The contents of <tt>yt</tt> are interpreted as in the function
 *  @ref forward_clusterbasis_avector :
 *  if <tt>cb</tt> is not a leaf, the first <tt>cb->koff[cb->directions]</tt> 
 *  rows of <tt>yt</tt> are coefficients to be multiplied by @f$V_{t \iota}@f$
 *  for all directions @f$ \iota \in \mathcal{D}_{t} @f$.
 *  Successively followed by <tt>cb->son[j]->k[iota1]</tt> coefficients with 
 *  @f$ iota1 = \mathrm{sd}_{t}(\iota) @f$ for all sons.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->koff[cb->directions]</tt> rows of
 *  <tt>yt</tt> are also coefficients to be multiplied by @f$V_{t \iota}@f$
 *  for all  directions @f$ \iota @f$,
 *  followed by <tt>cb->t->size</tt> coefficients in cluster numbering
 *  to be added to appropriate entries of the target vector.
 *
 *  @remark This function accesses the target vector <tt>y</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>y</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in
 *  the original numbering.
 *
 *  @param cb Directional cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector. */

HEADER_PREFIX void
backward_dclusterbasis(pcdclusterbasis cb, pavector yt, pavector y);

/** @brief Slow version of the forward transformation.
 *
 *  Version of @ref forward_dclusterbasis that uses the matrices
 *  <tt>cb->V</tt> for all clusters, not only for leaves.
 *  If <tt>cb</tt> is not a leaf cluster <tt>V</tt> is reconstructed.
 *
 *  @param cb Directional cluster basis.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */

HEADER_PREFIX void
slowforward_dclusterbasis(pcdclusterbasis cb, pcavector x, pavector xt);

/** @brief Slow version of the backward transformation.
 *
 *  Version of @ref backward_dclusterbasis that uses the matrices
 *  <tt>cb->V</tt> for all clusters, not only for leaves.
 *  If <tt>cb</tt> is not a leaf cluster <tt>V</tt> is reconstructed.
 *
 *  @param cb Directional cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector. */

HEADER_PREFIX void
slowbackward_dclusterbasis(pcdclusterbasis cb, pcavector yt, pavector y);

/* ------------------------------------------------------------
 * Forward and backward transformation for the root only
 * ------------------------------------------------------------ */

/** @brief Compute @f$\hat x_t = V_t^* x@f$.
 *
 *  This function computes @f$\hat x_t = V_t^* x@f$ using 
 *  the transfer matrices @f$E_{t \iota} @f$ and the leaf matrices @f$V_{t \iota}@f$,
 *  for all directions @f$ \iota @f$.
 *  Similar to @ref forward_dclusterbasis without permutation. 
 *
 *  @param cb Directional cluster basis.
 *  @param xp Source vector using cluster numbering corresponding
 *         to <tt>cb->t</tt>.
 *  @param xt Target vector of dimension <tt>cb->kbranch</tt>.
 *         Its first <tt>cb->koff[cb->directions]</tt>  rows will be filled by the
 *         transformed coefficients. */

HEADER_PREFIX void
compress_dclusterbasis(pcdclusterbasis cb,
		       pcavector xp, pavector xt);

/** @brief Add @f$ \alpha V_t \hat y_t@f$ to target vector @f$y@f$.
 *
 *  This function computes @f$y \gets y + V_t \hat y_t@f$ using 
 *  the transfer matrices @f$E_{t \iota} @f$ and the leaf matrices @f$V_{t \iota}@f$,
 *  for all directions @f$ \iota @f$.
 *  Similar to @ref backward_dclusterbasis without permutation. 
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param cb Directional cluster basis.
 *  @param yt Target vector of dimension <tt>cb->kbranch</tt>.
 *         Its first <tt>cb->koff[cb->directions]</tt> rows contain transformed coefficients.
 *  @param yp Source vector using cluster numbering corresponding
 *         to <tt>cb->t</tt>. */

HEADER_PREFIX void
expand_dclusterbasis(field alpha, pcdclusterbasis cb,
		     pavector yt, pavector yp);

/** @brief Compute @f$\widehat{X}_t = V_t^* X@f$.
 *
 *  Version of @ref compress_dclusterbasis that works for matrices
 *  instead of vectors, applying the transformation to all columns
 *  simultaneously.
 *
 *  @param cb Cluster basis.
 *  @param Xp Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows.
 *  @param Xt Target vector with <tt>cb->kbranch</tt> rows.
 *         Its first <tt>cb->k[cb->directions]</tt> rows will be filled by the
 *         transformed coefficients. */

HEADER_PREFIX void
blockcompress_dclusterbasis(pcdclusterbasis cb,
			    pcamatrix Xp, pamatrix Xt);

/** @brief Compute @f$ X =X + \alpha V_t \widehat{X}_t@f$.
 *
 *  Version of @ref expand_dclusterbasis that works for matrices
 *  instead of vectors, applying the transformation to all columns
 *  simultaneously.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param cb Cluster basis.
 *  @param Xt Source matrix using cluster numbering corresponding
 *         to <tt>cb->t</tt> in the rows.
 *  @param Xp Target vector with <tt>cb->kbranch</tt> rows.
 *         Its first <tt>cb->koff[cb->directions]</tt> rows will be filled by the
 *         transformed coefficients. */

HEADER_PREFIX void
blockexpand_dclusterbasis(field alpha, pcdclusterbasis cb,
			  pamatrix Xt, pamatrix Xp);

/* ------------------------------------------------------------
 * Enumeration by cluster number
 * ------------------------------------------------------------ */

/** @brief Enumerate @ref dclusterbasis according to 
 * the durectional cluster tree.
 *
 *  The @ref dclusterbasis elements are enumerated in an array of
 *  size <tt>t->desc</tt>. The enumeration starts with <tt>0</tt> assigned to
 *  the root and then proceeds with <tt>cb->son[0]</tt> corresponding to
 *  the entries <tt>1</tt> to <tt>t->son[0]->desc</tt> in the array and
 *  ending with <tt>cb->son[sons-1]</tt> corresponding to the last
 *  <tt>t->son[sons-1]->desc</tt> entries.
 *
 *  @param t Directional Cluster tree.
 *  @param cb Directional cluster basis matching the cluster tree given by <tt>t</tt>.
 *  @returns Array of size <tt>t->desc</tt> containing pointers to the
 *         @ref clusterbasis objects corresponding to <tt>cb</tt> and
 *         its descendants. */

HEADER_PREFIX pdclusterbasis *
enumerate_dclusterbasis(pcdcluster t, pdclusterbasis cb);

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

/** @brief Hierarchical iterator for a @ref dclusterbasis.
 *
 *  Iterate over the directional cluster basis and all its descendants.
 *  The <tt>pre</tt> function is called for each element before its
 *  descendants are processed, the <tt>post</tt> function is called
 *  afterwards.
 *
 *  <tt>cbname</tt> will take values between <tt>cbname</tt> and <tt>cbname+cb->t->sons-1,</tt>
 *  where <tt>cbname</tt> is interpreted as the index of the root <tt>cb,</tt>
 *  while the indices of the first son start at <tt>cbname+1,</tt> the indices
 *  of the second at <tt>cbname+1+cb->t->son[0]->desc,</tt> and so on.
 *
 *  @param cb Root of the directional cluster basis tree.
 *  @param tname Number of the root.
 *  @param pardepth Parallelization depth.
 *  @param pre Function to be called before the descendants of <tt>cb</tt>
 *         are processed.
 *  @param post Function to be called after the descentants of <tt>cb</tt>
 *         have been processed.
 *  @param data Additional data passed on to <tt>pre</tt> and <tt>post</tt>. */

HEADER_PREFIX void
iterate_dclusterbasis(pdclusterbasis cb, uint tname,
		      uint pardepth,
		      void (*pre)(pdclusterbasis cb, uint tname,
				  uint pardepth, void *data),
		      void (*post)(pdclusterbasis cb, uint tname,
				   uint pardepth, void *data),
		      void *data);

/* ------------------------------------------------------------
 * Orthogonalization
 * ------------------------------------------------------------ */


/** 
 * @brief Create an orthogonal directional cluster basis.
 *
 *  Compute an orthogonal directional cluster basis
 *  @f$(Q_t,\iota)_{t\in{\mathcal T}_{\mathcal I}, \iota\in D_t}@f$ 
 *  and factors @f$(R_t,\iota)_{t\in{\mathcal T}_{\mathcal I}\iota\in D_t}@f$
 *  such that @f$V_t,\iota = Q_t,\iota R_t,\iota@f$ holds for all directional
 *  clusters @f$t@f$ and directions @f$\iota \in D_{t}@f$.
 *
 *  @remark If the directional clusterbasis of an already existing directional
 *         @f$\mathcal{DH}^2@f$-matrix is orthogonalized, you need to call
 *         @ref resize_coupling_dh2matrix to get the new coupling matrices.
 *  @attention The @ref resize_coupling_dh2matrix needs the directional cluster operator
 *         filled in this function!
 *  @param cb Original directional cluster basis
 *         @f$(V_t,\iota)_{t\in{\mathcal T}_{\mathcal I}, \iota\in D_t}@f$,
 * 	        will be overwritten with the orthogonal cluster basis.
 *  @param co @ref dclusteroperator object will be overwritten with basis 
 *         change matrices.
 *         A simple way to set up this @ref dclusteroperator is
 *         to use @ref build_from_dclusterbasis_dclusteroperator.
 */

HEADER_PREFIX void
ortho_dclusterbasis(pdclusterbasis cb, pdclusteroperator co);

/** 
 * @brief Compute weight matrices with QR-decomposition.
 * 
 *  Compute factors @f$(R_t,\iota)_{t\in{\mathcal T}_{\mathcal I}\iota\in D_t}@f$
 *  for all directional clusters @f$t@f$ and save them in the directional
 *  cluster operator. 
 * 
 * @param cb Directional cluster basis.
 * @param co Directional cluster operator for saving the weight matrices.
 */

HEADER_PREFIX void
weight_dclusterbasis_dclusteroperator(pdclusterbasis cb, pdclusteroperator co);

/** 
 * @brief Check if a given directional cluster basis is orthogonal or not.
 *
 *  Evaluate for a directional leaf cluster
 *  @f$\|V_t,\iota^* V_t,\iota - I\|_F@f$ 
 *  and for a directional non-leaf cluster 
 *  @f$\|\widehat{V}_t,\iota^* \widehat{V}_t,\iota - I\|_F@f$
 *  and returns the maximum value.
 *  If the maximum is significantly larger than the machine eps,
 *  the directional cluster shouldn't be considered as orthogonal.
 *  
 *  @param cb Directional cluster basis which should be checked.
 *  @returns Returns a measure for the orthogonality of the @ref dclusterbasis
 *         <tt>cb</tt>. 
 */

HEADER_PREFIX real
check_ortho_dclusterbasis(pcdclusterbasis cb);


/** 
 * @brief Clone the structure of an already existing directional cluster basis.
 * 
 * Simply clone the structure without filling any matrices.
 * 
 * @param cb @ref dclusterbasis to be cloned.
 * @return Clone of the directional cluster basis <tt>cb</tt>.
 */

HEADER_PREFIX pdclusterbasis
clone_structure_dclusterbasis(pcdclusterbasis cb);

/** 
 * @brief Duplicate a directional cluster basis.
 * 
 * Simply copy the structure and fill all matrices.
 * 
 * @param cb @ref dclusterbasis to be duplicated.
 * @return Directional cluster basis with the copy of <tt>cb</tt>.
 */

HEADER_PREFIX pdclusterbasis
duplicate_dclusterbasis(pcdclusterbasis cb);

/** @} */

#endif
